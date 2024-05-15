import math
import numpy as np
import ngsolve as ngs

from netgen.webgui import Draw, register_draw_type, WebGLScene, encodeData


def updatePMinMax( pmat, pmima=None ):
    pmima_new = [ngs.Vector(pmat[:,i], copy=False).MinMax(ignore_inf=True) for i in range(3)]

    if pmima is not None:
        for i in range(3):
            minmax = ( min(pmima[i][0], pmima_new[i][0]), max(pmima[i][1], pmima_new[i][1]) )
            pmima_new[i] = minmax

    return pmima_new

def getMinMax( vals, fmin=None, fmax=None ):
    funcmin,funcmax = ngs.Vector(vals, copy=False).MinMax(True)
    if fmin is not None:
        funcmin = min(funcmin, fmin)
    if fmax is not None:
        funcmax = max(funcmax, fmax)
    return funcmin, funcmax


@register_draw_type(ngs.Mesh, ngs.CoefficientFunction, ngs.GridFunction, ngs.comp.GridFunctionCoefficientFunction)
def GetData(obj, args, kwargs):
        regions = {}

        cf = None
        if isinstance(obj, ngs.GridFunction):
            cf = obj
            mesh = obj.space.mesh
            if len(args) and isinstance(args[0], (ngs.Region, ngs.Mesh)):
                mesh = args[0]
        elif isinstance(obj, ngs.CoefficientFunction):
            cf = obj
            if len(args)<1:
                raise RuntimeError("Cannot draw CoefficientFunction without mesh")
            mesh = args[0]
        else:
            # assume, mesh is given
            mesh = obj

        # got one region as input: draw only on region and region.Boundaries()
        if isinstance(mesh, ngs.comp.Region):
            region = mesh
            mesh = region.mesh
            regions[region.VB()] = region
            while region.VB() != ngs.BBND:
                bnd_reg = region.Boundaries()
                regions[bnd_reg.VB()] = bnd_reg
                region = bnd_reg

        # draw exactly on given list of regions (no regions.Boundaries())
        elif isinstance(mesh, list):
            regions = {}
            list_ = mesh
            mesh = None
            for reg in list_:
                if isinstance(reg, ngs.comp.Region):
                    mesh = reg.mesh
                    regions[reg.VB()] = reg
                if isinstance(reg, ngs.comp.VorB):
                    regions[reg] = reg
            if mesh is None:
                raise RuntimeError("missing mesh/region")

        # draw on whole mesh
        elif isinstance(mesh, ngs.comp.Mesh):
            for vb in [ngs.VOL, ngs.BND, ngs.BBND]:
                regions[vb] = vb
        else:
            raise RuntimeError("invalid input mesh")

        # fill with empty regions
        for vb in [ngs.VOL, ngs.BND, ngs.BBND]:
            if vb not in regions:
                regions[vb] = mesh.Region(vb, None)

        if 'intpoints' not in kwargs:
            kwargs['intpoints']  = None
        d = BuildRenderData(mesh, cf, order=kwargs['order'], draw_surf=kwargs['draw_surf'], draw_vol=kwargs['draw_vol'], intpoints=kwargs['intpoints'], deformation=kwargs['deformation'], regions=regions, objects=kwargs['objects'], nodal_p1=kwargs['nodal_p1'], settings=kwargs['settings'])

        if isinstance(cf, ngs.GridFunction) and len(cf.vecs)>1:
            # multidim gridfunction - generate data for each component
            gf = ngs.GridFunction(cf.space)
            dim = len(cf.vecs)

            md_deformation = False;
            if 'deformation' in kwargs:
                deformation = kwargs['deformation']
                if isinstance(deformation, ngs.GridFunction) and len(deformation.vecs)==dim:
                    md_deformation = True
                    deformation = ngs.GridFunction(deformation.space)

            data = []
            for i in range(1,dim):
                gf.vec.data = cf.vecs[i]

                if md_deformation:
                    deformation.vec.data = deformation.vecs[i]

                data.append(BuildRenderData(mesh, gf, order=kwargs['order'], draw_surf=kwargs['draw_surf'], draw_vol=kwargs['draw_vol'], intpoints=kwargs['intpoints'], deformation=deformation, regions=regions, objects=kwargs['objects'], nodal_p1=kwargs['nodal_p1'], settings=kwargs['settings']))
            d['multidim_data'] = data
            d['multidim_interpolate'] = kwargs['interpolate_multidim']
            d['multidim_animate'] = kwargs['animate']

        d['deformation_scale'] = kwargs['scale']

        if 'is_complex' in d and d['is_complex'] and 'animate_complex' in kwargs:
            s = d['gui_settings']
            if 'Complex' not in s:
                s['Complex'] = dict(phase= 0.0, speed=2, animate=False)
            s['Complex']['animate'] = kwargs['animate_complex']

        if 'colors' in kwargs:
            d['colors'] = kwargs['colors']


        return d


bezier_trig_trafos = { }  # cache trafos for different orders

timer = ngs.Timer("BuildRenderData")
timerp1 = ngs.Timer("GetNodalP1Data")
timer2map = ngs.Timer("edges map")
timer2 = ngs.Timer("edges")
timermult = ngs.Timer("timer2 - mult")
timer3 = ngs.Timer("els")
timer3Bvals = ngs.Timer("timer3, bezier")
timer3minmax = ngs.Timer("els minmax")
timer2list = ngs.Timer("timer2 - make list")
timer3list = ngs.Timer("timer3 - make list")
timer4 = ngs.Timer("func")

timer3multnumpy = ngs.Timer("timer3 mul numpy")
timer3multngs = ngs.Timer("timer3 mul ngs")

def GetNodalP1Data(encode, mesh, cf, cf2=None):
    if cf2 is not None:
        cf = ngs.CF((cf,cf2))
    timerp1.Start()
    fes = ngs.NodalFESpace(mesh, order=1)**cf.dim
    gfu = ngs.GridFunction(fes)
    gfu.Interpolate(cf, ngs.BND)
    function_values = gfu.vec.FV().NumPy()
    nvert = mesh.nv
    function_values = function_values.reshape(-1, nvert).transpose().flatten()
    fmin, fmax = ngs.Vector(function_values, copy=False).MinMax(True)
    vertices = np.array(mesh.ngmesh._getVertices(), dtype=np.float32)

    pmat = vertices.reshape(-1, 3)
    pmin = pmat.min(axis=0)
    pmax = pmat.max(axis=0)
    mesh_center = list(np.array(0.5*(pmin+pmax), dtype=np.float64))
    mesh_radius = float(np.linalg.norm(pmax-pmin)/2)

    segments = mesh.ngmesh._getSegments()

    d = {}
    d['vertices'] = encode(vertices, dtype=np.float32)
    d['nodal_function_values'] = encode(function_values, dtype=np.float32)
    d['trigs'] = encode(np.array(mesh.ngmesh._get2dElementsAsTriangles(), dtype=np.int32))
    d['tets'] = encode(np.array(mesh.ngmesh._get3dElementsAsTets(), dtype=np.int32))
    d['segs'] = encode(np.array(segments, dtype=np.int32))
    d['mesh_center'] = mesh_center
    d['mesh_radius'] = mesh_radius
    d['funcmin'] = fmin
    d['funcmax'] = fmax

    timerp1.Stop()
    return d
    
    
def BuildRenderData(mesh, func, order=2, draw_surf=True, draw_vol=True, intpoints=None, deformation=None, regions=None, objects=[], nodal_p1=False, encoding='b64', settings={}):
    timer.Start()

    import inspect
    # stay backwards-compatible with older verisons of webgui_jupyter_widgets
    if 'encoding' in inspect.signature(encodeData).parameters:
        def encode(*args, **kwargs):
            return encodeData(*args, **kwargs, encoding=encoding)
    else:
        def encode(*args, **kwargs):
            return encodeData(*args, **kwargs)

    if isinstance(deformation, ngs.CoefficientFunction) and deformation.dim==2:
        deformation = ngs.CoefficientFunction((deformation, 0.0))

    #TODO: handle quads and non-smooth functions
    #TODO: subdivision

    d = {}
    d['gui_settings'] = settings
    d['ngsolve_version'] = ngs.__version__
    d['mesh_dim'] = mesh.dim
    # order = order or mesh.GetCurveOrder()
    if (not func) and (mesh.GetCurveOrder()==1) and (mesh.nv==len(mesh.ngmesh.Points())):
        order=1
    order2d = min(order, 3)
    order3d = min(order, 2)
    d['order2d'] = order2d
    d['order3d'] = order3d

    d['draw_vol'] = func and mesh.dim==3 and draw_vol and mesh.ne>0
    d['draw_surf'] = func and draw_surf

    d['objects'] = []
    for obj in objects:
        if isinstance(obj, dict):
            d['objects'].append(obj)
        else:
            d['objects'].append(obj._GetWebguiData())

    if isinstance(deformation, bool):
        d['deformation'] = deformation
        deformation = None

    func0 = None
    func2 = None
    if func and func.is_complex:
        d['is_complex'] = True
        func1 = func[0].real
        func2 = ngs.CoefficientFunction( (func[0].imag, 0.0) )
        d['funcdim'] = 2
    elif func and func.dim>1:
        func1 = func[0]
        func2 = ngs.CoefficientFunction( tuple(func[i] if i<func.dim else 0.0 for i in range(1,3)) ) # max 3-dimensional functions
        d['funcdim'] = func.dim
    elif func:
        func1 = func
        d['funcdim'] = 1
    else:
        # no function at all -> we are just drawing a mesh, eval mesh element index instead
        mats = mesh.GetMaterials()
        bnds = mesh.GetBoundaries()
        bbnds = mesh.GetBBoundaries()
        nmats = len(mesh.GetMaterials())
        nbnds = len(mesh.GetBoundaries())
        n = max(nmats, nbnds, len(bbnds))
        func1 = ngs.CoefficientFunction(list(range(n)))
        n_regions = [0, 0, nmats, nbnds]
        d['mesh_regions_2d'] = n_regions[mesh.dim]
        d['mesh_regions_3d'] = nmats if mesh.dim==3 else 0
        d['names'] = bnds if mesh.dim==3 else mats
        d['edge_names'] = bbnds if mesh.dim==3 else bnds
        d['funcdim'] = 0
        fds = mesh.ngmesh.FaceDescriptors()
        d["colors"] = [fd.color + (fd.transparency,) for fd in fds]

        func0 = func1

    if func0 is None:
        func0 = ngs.CoefficientFunction( 0.0 )

    d['show_wireframe'] = order2d>0
    d['show_mesh'] = order2d>0

    if order==1 and nodal_p1:
        d.update(GetNodalP1Data(encode, mesh, func1, func2))
        if '_override_data' in settings:
            d.update(settings['_override_data'])
        timer.Stop()
        return d

    func1 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, func1 ) )
    func0 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, func0 ) )

    if deformation is not None:
        func1 += ngs.CoefficientFunction((deformation, 0.0))
        func0 += ngs.CoefficientFunction((deformation, 0.0))

    if order2d>0:
        og = order2d
        timer2.Start()

        timer3Bvals.Start()

        # transform point-values to Bernsteinbasis
        def Binomial(n,i): return math.factorial(n) / math.factorial(i) / math.factorial(n-i)
        def Bernstein(x, i, n): return Binomial(n,i) * x**i*(1-x)**(n-i)
        Bvals = ngs.Matrix(og+1,og+1)
        for i in range(og+1):
            for j in range(og+1):
                Bvals[i,j] = Bernstein(i/og, j, og)
        iBvals = Bvals.I
        timer3Bvals.Stop()        
        # print (Bvals)
        # print (iBvals)

                
        Bezier_points = [] 

        # TODO: Quads
        ipts = [(i/og,0) for i in range(og+1)] + [(0, i/og) for i in range(og+1)] + [(i/og,1.0-i/og) for i in range(og+1)]
        ir_trig = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        ipts = [(i/og,0) for i in range(og+1)] + [(0, i/og) for i in range(og+1)] + [(i/og,1.0) for i in range(og+1)] + [(1.0, i/og) for i in range(og+1)]
        ir_quad = ngs.IntegrationRule(ipts, [0,]*len(ipts))

        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        reg = regions[vb]
        cf = func1 if draw_surf else func0
        timer2map.Start()
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, reg)
        timer2map.Stop()
        pmat = cf(pts)

        pmima = updatePMinMax(pmat)

        timermult.Start()
        pmat = pmat.reshape(-1, og+1, 4)
        if False:
            BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        else:
            BezierPnts = np.zeros( (og+1, pmat.shape[0], 4) )
            for i in range(4):
                ngsmat = ngs.Matrix(pmat[:,:,i].transpose()) 
                BezierPnts[:,:,i] = iBvals * ngsmat
        timermult.Stop()
        
        timer2list.Start()        
        for i in range(og+1):
            Bezier_points.append(encode(BezierPnts[i], dtype=np.float32))
        timer2list.Stop()        

        if func2 and draw_surf:
            pmat = func2(pts)
            pmat = pmat.reshape(-1, og+1, 2)
            timermult.Start()
            BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
            timermult.Stop()
            timer2list.Start()        
            for i in range(og+1):
                Bezier_points.append(encode(BezierPnts[i], dtype=np.float32))
            timer2list.Stop()        

        d['Bezier_points'] = Bezier_points

        ipts = [(i/og,0) for i in range(og+1)]
        ir_seg = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        vb = [ngs.VOL, ngs.BND, ngs.BBND][mesh.dim-1]
        reg = regions[vb]
        pts = mesh.MapToAllElements(ir_seg, reg)
        pmat = func0(pts)
        pmima = updatePMinMax(pmat)
        pmat = pmat.reshape(-1, og+1, 4)
        edge_data = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        edges = []
        for i in range(og+1):
            edges.append(encode(edge_data[i], dtype=np.float32))
        d['edges'] = edges

        timer2.Stop()
        timer3.Start()
        
        ndtrig = int((og+1)*(og+2)/2)
        
        if og in bezier_trig_trafos.keys():
            iBvals_trig = bezier_trig_trafos[og]
        else:
            def BernsteinTrig(x, y, i, j, n):
                return math.factorial(n)/math.factorial(i)/math.factorial(j)/math.factorial(n-i-j) \
                  * x**i*y**j*(1-x-y)**(n-i-j)
            Bvals = ngs.Matrix(ndtrig, ndtrig)
            ii = 0
            for ix in range(og+1):
                for iy in range(og+1-ix):
                    jj = 0
                    for jx in range(og+1):
                        for jy in range(og+1-jx):
                            Bvals[ii,jj] = BernsteinTrig(ix/og, iy/og, jx, jy, og)
                            jj += 1
                    ii += 1
            iBvals_trig = Bvals.I
            bezier_trig_trafos[og] = iBvals_trig

            
        # Bezier_points = [ [] for i in range(ndtrig) ]
        Bezier_points = []
        
        ipts = [(i/og,j/og) for j in range(og+1) for i in range(og+1-j)]
        ir_trig = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        ipts = ([(i/og,j/og) for j in range(og+1) for i in range(og+1-j)] + 
                [(1-i/og,1-j/og) for j in range(og+1) for i in range(og+1-j)])
        ir_quad = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        
        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        reg = regions[vb]
        if intpoints is not None:
            pts = mesh.MapToAllElements(intpoints, reg)
        else:
            pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, reg)

        pmat = ngs.CoefficientFunction( func1 if draw_surf else func0 ) (pts)

        
        timer3minmax.Start()
        pmima = updatePMinMax(pmat, pmima)
        funcmin,funcmax = getMinMax(pmat[:,3])
        timer3minmax.Stop()

        pmin, pmax = [ngs.Vector(p) for p in zip(*pmima)]
        mesh_center = 0.5*(pmin+pmax)
        mesh_radius = np.linalg.norm(pmax-pmin)/2
        
        pmat = pmat.reshape(-1, len(ir_trig), 4)

        if False:
            timer3multnumpy.Start()
            BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))
            timer3multnumpy.Stop()
        else:
            timer3multngs.Start()
            BezierPnts = np.zeros( (len(ir_trig), pmat.shape[0], 4) )
            for i in range(4):
                ngsmat = ngs.Matrix(pmat[:,:,i].transpose())
                BezierPnts[:,:,i] = iBvals_trig * ngsmat
            timer3multngs.Stop()

        timer3list.Start()        
        for i in range(ndtrig):
            Bezier_points.append(encode(BezierPnts[i], dtype=np.float32))
        timer3list.Stop()        

        if func2 and draw_surf:
            pmat = ngs.CoefficientFunction( func2 ) (pts)   

            pmat = pmat.reshape(-1, len(ir_trig), 2)

            funcmin, funcmax = getMinMax(pmat.flatten(), funcmin, funcmax)
            BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))
            if og==1:
                for i in range(ndtrig):
                    Bezier_points.append(encode(BezierPnts[i], dtype=np.float32))
            else:
                BezierPnts = BezierPnts.transpose((1,0,2)).reshape(-1, len(ir_trig)//2, 4).transpose((1,0,2))

                for i in range(ndtrig//2):
                    Bezier_points.append(encode(BezierPnts[i], dtype=np.float32))

        d['Bezier_trig_points'] = Bezier_points    
        d['mesh_center'] = list(mesh_center)
        d['mesh_radius'] = mesh_radius
        timer3.Stop()



    timer4.Start()

    if d['draw_vol']:
        p0 = []
        p1 = []
        p2 = []
        p3 = []
        values = []
        tets = []

        midpoint = lambda p0, p1: tuple((0.5*(p0[i]+p1[i]) for i in range(3)))

        def makeP2Tets( p1_tets ):
            p2_tets = []
            for tet in p1_tets:
                tet.append( midpoint(tet[0], tet[3]) )
                tet.append( midpoint(tet[1], tet[3]) )
                tet.append( midpoint(tet[2], tet[3]) )
                tet.append( midpoint(tet[0], tet[1]) )
                tet.append( midpoint(tet[0], tet[2]) )
                tet.append( midpoint(tet[1], tet[2]) )
                p2_tets.append(tet)
            return p2_tets

        # divide any element into tets
        p1_tets = {}
        p1_tets[ngs.ET.TET]   = [[(1,0,0), (0,1,0), (0,0,1), (0,0,0)]]
        p1_tets[ngs.ET.PYRAMID]=[[(1,0,0), (0,1,0), (0,0,1), (0,0,0)],
                                 [(1,0,0), (0,1,0), (0,0,1), (1,1,0)]]
        p1_tets[ngs.ET.PRISM] = [[(1,0,0), (0,1,0), (0,0,1), (0,0,0)],
                                 [(0,0,1), (0,1,0), (0,1,1), (1,0,0)],
                                 [(1,0,1), (0,1,1), (1,0,0), (0,0,1)]]
        p1_tets[ngs.ET.HEX]   = [[(1,0,0), (0,1,0), (0,0,1), (0,0,0)],
                                 [(0,1,1), (1,1,1), (1,1,0), (1,0,1)],
                                 [(1,0,1), (0,1,1), (1,0,0), (0,0,1)],
                                 [(0,1,1), (1,1,0), (0,1,0), (1,0,0)],
                                 [(0,0,1), (0,1,0), (0,1,1), (1,0,0)],
                                 [(1,0,1), (1,1,0), (0,1,1), (1,0,0)]]

        intrules = {}
        for eltype in p1_tets:
          points = p1_tets[eltype]
          if order3d>1:
            points = makeP2Tets( points )
          intrules[eltype] = ngs.IntegrationRule( sum(points, []) )

        if intpoints is not None:
            pts = mesh.MapToAllElements(intpoints, regions[ngs.VOL])
        else:
            pts = mesh.MapToAllElements(intrules, regions[ngs.VOL])
            
        pmat = func1(pts)

        np_per_tet = len(intrules[ngs.ET.TET])

        ne = mesh.GetNE(ngs.VOL)
        pmat = pmat.reshape(-1, np_per_tet, 4)
        
        funcmin, funcmax = getMinMax(pmat[:,:,3].flatten(), funcmin, funcmax)
        points3d = []
        for i in range(np_per_tet):
            points3d.append(encode(pmat[:,i,:], dtype=np.float32))

        if func2:
            pmat = func2(pts).reshape(-1, np_per_tet//2, 4)
            funcmin, funcmax = getMinMax(pmat.flatten(), funcmin, funcmax)
            for i in range(np_per_tet//2):
                points3d.append(encode(pmat[:,i,:], dtype=np.float32))
        d['points3d'] = points3d
    if func:
        d['funcmin'] = funcmin
        d['funcmax'] = funcmax
    if '_override_data' in settings:
        d.update(settings['_override_data'])
    timer4.Stop()
    timer.Stop()
    return d

__all__ = ['Draw']
