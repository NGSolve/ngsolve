import math
import numpy as np
from time import time
import ngsolve as ngs
import os

from webgui_jupyter_widgets import BaseWebGuiScene, encodeData, WebGuiDocuWidget
import webgui_jupyter_widgets.widget as wg

class WebGLScene(BaseWebGuiScene):
    def __init__(self, cf, mesh, order, min_, max_, draw_vol, draw_surf, autoscale, deformation, interpolate_multidim, animate, clipping, vectors, on_init, eval_function, eval_):
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

        self.deformation = deformation

        if isinstance(mesh, ngs.comp.Region):
            self.region = mesh
            self.mesh = self.region.mesh
        else:
            self.region = None

    def GetData(self, set_minmax=True):
        import json
        d = BuildRenderData(self.mesh, self.cf, self.order, draw_surf=self.draw_surf, draw_vol=self.draw_vol, deformation=self.deformation, region=self.region)

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

                data.append(BuildRenderData(self.mesh, gf, self.order, draw_surf=self.draw_surf, draw_vol=self.draw_vol, deformation=deformation, region=self.region))
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
timer2 = ngs.Timer("edges")
timermult = ngs.Timer("timer2 - mult")
timer3 = ngs.Timer("els")
timer3Bvals = ngs.Timer("timer3, bezier")
timer3minmax = ngs.Timer("els minmax")
timer2list = ngs.Timer("timer2 - make list")
timer3list = ngs.Timer("timer3 - make list")
timer4 = ngs.Timer("func")

    
def BuildRenderData(mesh, func, order=2, draw_surf=True, draw_vol=True, deformation=None, region=True):
    timer.Start()

    if isinstance(deformation, ngs.CoefficientFunction) and deformation.dim==2:
        deformation = ngs.CoefficientFunction((deformation, 0.0))

    #TODO: handle quads and non-smooth functions
    #TODO: subdivision

    d = {}
    d['ngsolve_version'] = ngs.__version__
    d['mesh_dim'] = mesh.dim
    # order = order or mesh.GetCurveOrder()
    if (not func) and (mesh.GetCurveOrder()==1):
        order=1
    order2d = min(order, 3)
    order3d = min(order, 2)
    d['order2d'] = order2d
    d['order3d'] = order3d

    d['draw_vol'] = func and mesh.dim==3 and draw_vol and mesh.ne>0
    d['draw_surf'] = func and draw_surf

    if isinstance(deformation, bool):
        d['deformation'] = deformation
        deformation = None

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
        nmats = len(mesh.GetMaterials())
        nbnds = len(mesh.GetBoundaries())
        n = max(nmats, nbnds)
        func1 = ngs.CoefficientFunction(list(range(n)))
        n_regions = [0, 0, nmats, nbnds]
        d['mesh_regions_2d'] = n_regions[mesh.dim]
        d['mesh_regions_3d'] = nmats if mesh.dim==3 else 0
        d['funcdim'] = 0
    func1 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, func1 ) )
    func0 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, 0.0 ) )
    if deformation is not None:
        func1 += ngs.CoefficientFunction((deformation, 0.0))
        func0 += ngs.CoefficientFunction((deformation, 0.0))

    d['show_wireframe'] = False
    d['show_mesh'] = False
    if order2d>0:
        og = order2d
        d['show_wireframe'] = True
        d['show_mesh'] = True
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
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)
        pmat = cf(pts)

        timermult.Start()
        pmat = pmat.reshape(-1, og+1, 4)
        BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        timermult.Stop()
        
        timer2list.Start()        
        for i in range(og+1):
            Bezier_points.append(encodeData(BezierPnts[i]))
        timer2list.Stop()        

        if func2 and draw_surf:
            pmat = func2(pts)
            pmat = pmat.reshape(-1, og+1, 2)
            timermult.Start()
            BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
            timermult.Stop()
            timer2list.Start()        
            for i in range(og+1):
                Bezier_points.append(encodeData(BezierPnts[i]))
            timer2list.Stop()        

        d['Bezier_points'] = Bezier_points

        ipts = [(i/og,0) for i in range(og+1)]
        ir_seg = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        vb = [ngs.VOL, ngs.BND, ngs.BBND][mesh.dim-1]
        if region and region.VB() == vb:
            vb = region
        pts = mesh.MapToAllElements(ir_seg, vb)
        pmat = func0(pts)
        pmat = pmat.reshape(-1, og+1, 4)
        edge_data = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        edges = []
        for i in range(og+1):
            edges.append(encodeData(edge_data[i]))
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
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)

        pmat = ngs.CoefficientFunction( func1 if draw_surf else func0 ) (pts)
        timer3minmax.Start()
        funcmin = np.min(pmat[:,3])
        funcmax = np.max(pmat[:,3])
        pmin = np.min(pmat[:,0:3], axis=0)
        pmax = np.max(pmat[:,0:3], axis=0)
        mesh_center = (pmin+pmax)/2
        mesh_radius = np.linalg.norm(pmax-pmin)/2
        timer3minmax.Stop()

        pmat = pmat.reshape(-1, len(ir_trig), 4)
        BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))

        timer3list.Start()        
        for i in range(ndtrig):
            Bezier_points.append(encodeData(BezierPnts[i]))
        timer3list.Stop()        

        if func2 and draw_surf:
            pmat = ngs.CoefficientFunction( func2 ) (pts)   

            pmat = pmat.reshape(-1, len(ir_trig), 2)

            funcmin = min(funcmin, np.min(pmat))
            funcmax = max(funcmax, np.max(pmat))
            BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))
            if og==1:
                for i in range(ndtrig):
                    Bezier_points.append(encodeData(BezierPnts[i]))
            else:
                BezierPnts = BezierPnts.transpose((1,0,2)).reshape(-1, len(ir_trig)//2, 4).transpose((1,0,2))

                for i in range(ndtrig//2):
                    Bezier_points.append(encodeData(BezierPnts[i]))

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

        if order3d==1:
            ipts = [(1,0,0), (0,1,0), (0,0,1), (0,0,0)]
            ir_tet = ngs.IntegrationRule( ipts, [0]*len(ipts) )


            ipts = ([(1,0,0), (0,1,0), (0,0,1), (0,0,0)] +
                         [(0,0,1), (0,1,0), (0,1,1), (1,0,0)] +
                         [(1,0,1), (0,1,1), (1,0,0), (0,0,1)])
            ir_prism = ngs.IntegrationRule( ipts, [0]*len(ipts) )


            ipts_hex = ([(1,0,0), (0,1,0), (0,0,1), (0,0,0)] +
                         [(0,1,1), (1,1,1), (1,1,0), (1,0,1)] +
                         [(1,0,1), (0,1,1), (1,0,0), (0,0,1)] +
                         [(0,1,1), (1,1,0), (0,1,0), (1,0,0)] +
                         [(0,0,1), (0,1,0), (0,1,1), (1,0,0)] +
                         [(1,0,1), (1,1,0), (0,1,1), (1,0,0)] )
            ir_hex = ngs.IntegrationRule( ipts_hex, [0]*len(ipts_hex) )
            pts = mesh.MapToAllElements({ngs.ET.TET: ir_tet, ngs.ET.PRISM: ir_prism, ngs.ET.HEX: ir_hex }, ngs.VOL)
            
            
        else:
            # TODO: make prism
            ir_tet = ngs.IntegrationRule( [
                (1,0,0),
                (0,1,0),
                (0,0,1),
                (0,0,0),
                (0.5,0,0),
                (0,0.5,0),
                (0,0,0.5),
                (0.5,0.5,0),
                (0.5,0,0.5),
                (0,0.5,0.5) ],
                [0]*10 )
            pts = mesh.MapToAllElements({ngs.ET.TET: ir_tet}, ngs.VOL)

        pmat = func1(pts)

        ne = mesh.GetNE(ngs.VOL)
        pmat = pmat.reshape(-1, len(ir_tet), 4)
        
        funcmin = min(funcmin, np.min(pmat[:,:,3]))
        funcmax = max(funcmax, np.max(pmat[:,:,3]))
        points3d = []
        for i in range(len(ir_tet)):
            points3d.append(encodeData(pmat[:,i,:]))

        if func2:
            pmat = func2(pts).reshape(-1, len(ir_tet)//2, 4)
            funcmin = min(funcmin, np.min(pmat))
            funcmax = max(funcmax, np.max(pmat))
            for i in range(len(ir_tet)//2):
                points3d.append(encodeData(pmat[:,i,:]))
        d['points3d'] = points3d
    if func:
        d['funcmin'] = funcmin
        d['funcmax'] = funcmax
    timer4.Stop()
    timer.Stop()
    return d

def Draw(mesh_or_func, mesh_or_none=None, name='function', order=2, min=None, max=None, draw_vol=True, draw_surf=True, autoscale=True, deformation=False, interpolate_multidim=False, animate=False, clipping=None, vectors=None, js_code=None, eval_function=None, eval=None, filename=""):
    if isinstance(mesh_or_func, ngs.Mesh):
        mesh = mesh_or_func
        func = None

    if isinstance(mesh_or_func, ngs.CoefficientFunction):
        func = mesh_or_func
        mesh = mesh_or_none

    if isinstance(mesh_or_func, ngs.GridFunction):
        func = mesh_or_func
        mesh = mesh_or_none or func.space.mesh
        
    scene = WebGLScene(func, mesh, order, min_=min, max_=max, draw_vol=draw_vol, draw_surf=draw_surf, autoscale=autoscale, deformation=deformation, interpolate_multidim=interpolate_multidim, animate=animate, clipping=clipping, vectors=vectors, on_init=js_code, eval_function=eval_function, eval_=eval)
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


def _DrawDocu(mesh_or_func, mesh_or_none=None, name='function', order=2, min=None, max=None, draw_vol=True, draw_surf=True, autoscale=True, deformation=False, interpolate_multidim=False, animate=False, clipping=None, vectors=None, js_code=None, eval_function=None, eval=None, filename=""):
    if isinstance(mesh_or_func, ngs.Mesh):
        mesh = mesh_or_func
        func = None

    if isinstance(mesh_or_func, ngs.CoefficientFunction):
        func = mesh_or_func
        mesh = mesh_or_none

    if isinstance(mesh_or_func, ngs.GridFunction):
        func = mesh_or_func
        mesh = mesh_or_none or func.space.mesh
        
    scene = WebGLScene(func, mesh, order, min_=min, max_=max, draw_vol=draw_vol, draw_surf=draw_surf, autoscale=autoscale, deformation=deformation, interpolate_multidim=interpolate_multidim, animate=animate, clipping=clipping, vectors=vectors, on_init=js_code, eval_function=eval_function, eval_=eval)
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

