import math
import numpy as np
from time import time
import ngsolve as ngs
import os

from webgui_jupyter_widgets import BaseWebGuiScene, encodeData, WebGuiDocuWidget
import webgui_jupyter_widgets.widget as wg

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


class WebGLScene(BaseWebGuiScene):
    def __init__(self, cf, mesh, order, min_, max_, draw_vol, draw_surf, autoscale, intpoints, deformation, interpolate_multidim, animate, clipping, vectors, on_init, eval_function, eval_, objects, nodal_p1, settings={}):
        from IPython.display import display, Javascript
        import threading
        self.cf = cf
        self.mesh = mesh
        self.order = order
        self.min = min_
        self.max = max_
        self.draw_vol = draw_vol
        self.draw_surf = draw_surf
        self.autoscale = autoscale
        self.interpolate_multidim = interpolate_multidim
        self.animate = animate
        self.clipping = clipping
        self.vectors = vectors
        self.on_init = on_init
        self.eval_function = eval_function
        self.eval_ = eval_
        self.objects = objects
        self.nodal_p1 = nodal_p1
        self.settings = settings

        self.intpoints = intpoints
        self.deformation = deformation
        self.encoding = 'b64'

        if isinstance(mesh, ngs.comp.Region):
            self.region = mesh
            self.mesh = self.region.mesh
        else:
            self.region = None

    def GetData(self, set_minmax=True):
        encoding = self.encoding
        import json
        d = BuildRenderData(self.mesh, self.cf, self.order, draw_surf=self.draw_surf, draw_vol=self.draw_vol, intpoints=self.intpoints, deformation=self.deformation, region=self.region, objects=self.objects, nodal_p1=self.nodal_p1, encoding=encoding, settings=self.settings)

        if isinstance(self.cf, ngs.GridFunction) and len(self.cf.vecs)>1:
            # multidim gridfunction - generate data for each component
            gf = ngs.GridFunction(self.cf.space)
            dim = len(self.cf.vecs)

            if isinstance(self.deformation, ngs.GridFunction) and len(self.deformation.vecs)==dim:
                md_deformation = True
                deformation = ngs.GridFunction(self.deformation.space)
            else:
                md_deformation = False
                deformation = self.deformation

            data = []
            for i in range(1,dim):
                gf.vec.data = self.cf.vecs[i]

                if md_deformation:
                    deformation.vec.data = self.deformation.vecs[i]

                data.append(BuildRenderData(self.mesh, gf, self.order, draw_surf=self.draw_surf, draw_vol=self.draw_vol, intpoints=self.intpoints, deformation=deformation, region=self.region, objects=self.objects, nodal_p1=self.nodal_p1, encoding=encoding, settings=self.settings))
            d['multidim_data'] = data
            d['multidim_interpolate'] = self.interpolate_multidim
            d['multidim_animate'] = self.animate


        if set_minmax:
            if self.min is not None:
                d['funcmin'] = self.min
            if self.max is not None:
                d['funcmax'] = self.max
            d['autoscale'] = self.autoscale

        if self.clipping is not None:
            d['clipping'] = True
            if isinstance(self.clipping, dict):
                allowed_args = ("x", "y", "z", "dist", "function", "pnt", "vec")
                if "vec" in self.clipping:
                    vec = self.clipping["vec"]
                    self.clipping["x"] = vec[0]
                    self.clipping["y"] = vec[1]
                    self.clipping["z"] = vec[2]
                if "pnt" in self.clipping:
                    d['mesh_center'] = list(self.clipping["pnt"])
                for name, val in self.clipping.items():
                    if not (name in allowed_args):
                        raise Exception('Only {} allowed as arguments for clipping!'.format(", ".join(allowed_args)))
                    d['clipping_' + name] = val

        if self.vectors is not None:
            d['vectors'] = True
            if isinstance(self.vectors, dict):
                for name, val in self.vectors.items():
                    if not (name in ("grid_size", "offset")):
                        raise Exception('Only "grid_size" and "offset" allowed as arguments for vectors!')
                    d['vectors_' + name] = val

        if self.on_init:
            d['on_init'] = self.on_init

        if self.eval_function:
            d['user_eval_function'] = self.eval_function

        # see shaders/utils.h for value explanation (function_mode)
        eval_ = self.eval_
        if eval_ is not None:
            if isinstance(eval_, int):
                d['eval'] = eval_
            elif eval_ == 'norm':
                d['eval'] = 3
            elif eval_ == 'real':
                d['eval'] = 5
            elif eval_ == 'imag':
                d['eval'] = 6

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
    
    
def BuildRenderData(mesh, func, order=2, draw_surf=True, draw_vol=True, intpoints=None, deformation=None, region=True, objects=[], nodal_p1=False, encoding='b64', settings={}):
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
        if region and region.VB() == vb:
            vb = region
        cf = func1 if draw_surf else func0
        timer2map.Start()
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)
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
        if region and region.VB() == vb:
            vb = region
        pts = mesh.MapToAllElements(ir_seg, vb)
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
        if region and region.VB() == vb:
            vb = region
        if intpoints is not None:
            pts = mesh.MapToAllElements(intpoints, vb)
        else:
            pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)

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
            pts = mesh.MapToAllElements(intpoints, ngs.VOL)            
        else:
            pts = mesh.MapToAllElements(intrules, ngs.VOL)
            
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

def Draw(mesh_or_func, mesh_or_none=None, name='function', order=2, min=None, max=None, draw_vol=True, draw_surf=True, autoscale=True, intpoints=None, deformation=False, interpolate_multidim=False, animate=False, clipping=None, vectors=None, js_code=None, eval_function=None, eval=None, filename="", objects=[], nodal_p1=False, settings={}):
    if isinstance(mesh_or_func, ngs.Mesh):
        mesh = mesh_or_func
        func = None

    if isinstance(mesh_or_func, ngs.CoefficientFunction):
        func = mesh_or_func
        mesh = mesh_or_none

    if isinstance(mesh_or_func, ngs.GridFunction):
        func = mesh_or_func
        mesh = mesh_or_none or func.space.mesh
        
    scene = WebGLScene(func, mesh, order, min_=min, max_=max, draw_vol=draw_vol, draw_surf=draw_surf, autoscale=autoscale, intpoints=intpoints, deformation=deformation, interpolate_multidim=interpolate_multidim, animate=animate, clipping=clipping, vectors=vectors, on_init=js_code, eval_function=eval_function, eval_=eval, objects=objects, nodal_p1=nodal_p1, settings=settings)
    if wg._IN_IPYTHON:
        if wg._IN_GOOGLE_COLAB:
            from IPython.display import display, HTML
            html = scene.GenerateHTML()
            display(HTML(html))
        else:
            # render scene using widgets.DOMWidget
            scene.Draw()
            return scene
    else:
        if filename:
            scene.GenerateHTML(filename=filename)
        return scene


def _DrawDocu(mesh_or_func, mesh_or_none=None, name='function', order=2, min=None, max=None, draw_vol=True, draw_surf=True, autoscale=True, intpoints=None, deformation=False, interpolate_multidim=False, animate=False, clipping=None, vectors=None, js_code=None, eval_function=None, eval=None, filename="", objects=[], nodal_p1=False, settings={}):
    if isinstance(mesh_or_func, ngs.Mesh):
        mesh = mesh_or_func
        func = None

    if isinstance(mesh_or_func, ngs.CoefficientFunction):
        func = mesh_or_func
        mesh = mesh_or_none

    if isinstance(mesh_or_func, ngs.GridFunction):
        func = mesh_or_func
        mesh = mesh_or_none or func.space.mesh
        
    scene = WebGLScene(func, mesh, order, min_=min, max_=max, draw_vol=draw_vol, draw_surf=draw_surf, autoscale=autoscale, intpoints=intpoints, deformation=deformation, interpolate_multidim=interpolate_multidim, animate=animate, clipping=clipping, vectors=vectors, on_init=js_code, eval_function=eval_function, eval_=eval, objects=objects, nodal_p1=nodal_p1, settings=settings)
    import json

    docu_path = os.environ['NETGEN_DOCUMENTATION_OUT_DIR']
    src_path = os.environ['NETGEN_DOCUMENTATION_SRC_DIR']
    cwd_path = os.path.abspath('.')
    rel_path = os.path.relpath('.', src_path)
    path = os.path.join(docu_path, rel_path)

    if not os.path.exists(path):
        os.makedirs(path)
    counter_file = os.path.join(docu_path, '.counter')
    if os.path.exists(counter_file):
        file_counter = int(open(counter_file,'r').read())+1
    else:
        file_counter = 0

    open(counter_file,'w').write(str(file_counter))

    data_file = 'render_data_{}.json'.format(file_counter)
    data_file_abs = os.path.join(path, data_file)
    preview_file = 'preview_{}.png'.format(file_counter)
    preview_file_abs = os.path.join(path, preview_file)


    widget = WebGuiDocuWidget()
    widget.value = {'render_data' : data_file, 'preview' : preview_file }
    scene.widget = widget
    data = scene.GetData()
    json.dump(data, open(data_file_abs, "w"))
    scene.MakeScreenshot(preview_file_abs, 1200, 600)
    scene.Redraw = lambda : None
    from IPython.display import display, HTML
    display(widget)
    return scene

if 'NETGEN_DOCUMENTATION_SRC_DIR' in os.environ:
    # we are buiding the documentation, some things are handled differently:
    # 1) Draw() is generating a .png (using headless chromium via selenium) and a render_data.json
    #    to show a preview image and load the render_data only when requested by user
    # 2) return a NGSDocuWebGuiWidget instead of NGSWebGuiWidget implementing the preview/load on demand of webgui

    _Draw = Draw
    Draw = _DrawDocu

