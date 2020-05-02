import math
import numpy as np
from time import time
import ngsolve as ngs

# the build script fills the contents of the variables below
shader_codes = {}
render_js_code = ""

websocket_port = 8765
html_template = """
<!DOCTYPE html>
<html>
    <head>
        <title>NGSolve WebGUI</title>
        <meta name='viewport' content='width=device-width, user-scalable=no'/>
        <style>
            body{
                margin:0;
                overflow:hidden;
            }
            canvas{
                cursor:grab;
                cursor:-webkit-grab;
                cursor:-moz-grab;
            }
            canvas:active{
                cursor:grabbing;
                cursor:-webkit-grabbing;
                cursor:-moz-grabbing;
            }
        </style>
    </head>
    <body>
          <script src="https://cdn.jsdelivr.net/npm/three@0.115.0/build/three.min.js"></script>
          <script src="https://cdn.jsdelivr.net/npm/three@0.115.0/examples/js/controls/OrbitControls.js"></script>
          <script src="https://cdnjs.cloudflare.com/ajax/libs/dat-gui/0.7.6/dat.gui.js"></script>
          <script src="https://cdnjs.cloudflare.com/ajax/libs/stats.js/r16/Stats.min.js"></script>
          {shader}
          <script>
            var render_data = {data}
            {render}
          </script>
    </body>
</html>
"""

# Copied from Ipython.display.IFrame, but using srcdoc instead of src
#   -> avoids the need to write an extra html-file
#   -> the html is automatically saved with the notebook and available in html exports
class MyIFrame(object):
    iframe = """
        <iframe
            width="{width}"
            height="{height}"
            srcdoc="{srcdoc}"
            frameborder="0"
            allowfullscreen
        ></iframe>
        """

    def __init__(self, srcdoc, width, height):
        self.srcdoc = srcdoc
        self.width = width
        self.height = height

    def _repr_html_(self):
        return self.iframe.format(srcdoc=self.srcdoc,
                                  width=self.width,
                                  height=self.height)

def readShaders():
    from glob import glob 
    import os.path
    codes = ""
    for name, code in shader_codes.items():
        if name.endswith(".h"):
            typ = "x-header"
        if name.endswith(".vert"):
            typ = "x-vertex"
        if name.endswith(".frag"):
            typ = "x-fragment"
        codes += '<script id={} type="x-shader/{}">\n{}</script>\n'.format(name, typ, code)
    return codes

def writeHTML(d):
    import json
    data = json.dumps(d)

    html = html_template.replace('{shader}', readShaders() )
    html = html.replace('{data}', data )
    html = html.replace('{render}', render_js_code )

    from IPython.display import display, IFrame

    ## write to html and use iframe<src="./output.html">
    open('output.html','w').write( html )
    # frame = IFrame(src='./output.html', width=700, height=400)

    # use <iframe srcdoc=html> to display html inline
    html = html.replace('"', '&quot;')
    frame = MyIFrame(srcdoc=html, width="100%", height=400)

    display(frame)
    # return frame

class WebGLScene:
    def __init__(self, cf, mesh, order, start_websocket=False):
        import threading
        self.cf = cf
        self.mesh = mesh
        self.order = order
        self.have_websocket = start_websocket
        if start_websocket:
            self.connected = set()
            self.thread = threading.Thread(target = self._eventLoopFunction)
            self.thread.start()

    def _send(self, data):
        import asyncio
        for websocket in self.connected.copy():
            asyncio.ensure_future(websocket.send(data))

    def Redraw(self):
        import numpy, json
        if self.have_websocket and len(self.connected):
            d = BuildRenderData(self.mesh, self.cf, self.order)
            data = json.dumps(d["Bezier_trig_points"])
            self.loop.call_soon_threadsafe(self._send, data.encode('ascii'))

    async def _handler(self, websocket, path):
        import websockets
        self.connected.add(websocket)
        try:
            while True:
                s = await websocket.recv()
                print("received websocket message", s)
        except websockets.exceptions.ConnectionClosed:
            pass
        finally:
            self.connected.remove(websocket)

    def _eventLoopFunction(self):
        import websockets, asyncio
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        self.server = websockets.serve(self._handler, "localhost", websocket_port)
        self.loop.run_until_complete(self.server)
        self.loop.run_forever()

    def __repr__(self):
        return ""


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

    
def BuildRenderData(mesh, func, order=2):

    
    timer.Start()

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
        func1 = ngs.CoefficientFunction(0.0)
        d['funcdim'] = 0
    func1 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, func1 ) )

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
        ir = ngs.IntegrationRule(ipts, [0,]*len(ipts))

        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        pts = mesh.MapToAllElements(ir, vb)
        pmat = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z) ) (pts)

        timermult.Start()
        pmat = pmat.reshape(-1, og+1, 3)
        BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        timermult.Stop()
        
        timer2list.Start()        
        for i in range(og+1):
            Bezier_points.append(encodeData(BezierPnts[i]))
        timer2list.Stop()        


            
        d['Bezier_points'] = Bezier_points

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
        
        ir = ngs.IntegrationRule( [(i/og,j/og) for j in range(og+1) for i in range(og+1-j)], [0,]*(ndtrig) )

        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        pts = mesh.MapToAllElements(ir, vb)
        pmat = ngs.CoefficientFunction( func1 ) (pts)
        timer3minmax.Start()
        funcmin = np.min(pmat[:,3])
        funcmax = np.max(pmat[:,3])
        pmin = np.min(pmat[:,0:3], axis=0)
        pmax = np.max(pmat[:,0:3], axis=0)
        mesh_center = (pmin+pmax)/2
        mesh_radius = np.linalg.norm(pmax-pmin)/2
        timer3minmax.Stop()
        pmat = pmat.reshape(mesh.GetNE(vb), len(ir), 4)
        BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))

        timer3list.Start()        
        for i in range(ndtrig):
            Bezier_points.append(encodeData(BezierPnts[i]))
        timer3list.Stop()        

        if func2:
            pmat = ngs.CoefficientFunction( func2 ) (pts)
            pmat = pmat.reshape(mesh.GetNE(vb), len(ir), 2)
            funcmin = min(funcmin, np.min(pmat))
            funcmax = max(funcmax, np.max(pmat))
            BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))
            if og==1:
                for i in range(ndtrig):
                    Bezier_points.append(encodeData(BezierPnts[i]))
            else:
                BezierPnts = BezierPnts.transpose((1,0,2)).reshape(mesh.GetNE(vb), len(ir)//2, 4).transpose((1,0,2))
                for i in range(ndtrig//2):
                    Bezier_points.append(encodeData(BezierPnts[i]))

        d['Bezier_trig_points'] = Bezier_points    
        d['mesh_center'] = list(mesh_center)
        d['mesh_radius'] = mesh_radius
        timer3.Stop()



    timer4.Start()

    if mesh.dim==3:
        p0 = []
        p1 = []
        p2 = []
        p3 = []
        values = []
        tets = []

        if order3d==1:
            ir = ngs.IntegrationRule( [(1,0,0), (0,1,0), (0,0,1), (0,0,0)], [0]*4 )
        else:
            ir = ngs.IntegrationRule( [
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
        pts = mesh.MapToAllElements(ir, ngs.VOL)
        pmat = func1(pts)

        ne = mesh.GetNE(ngs.VOL)
        pmat = pmat.reshape(ne, len(ir), 4)
        funcmin = min(funcmin, np.min(pmat[:,:,3]))
        funcmax = max(funcmax, np.max(pmat[:,:,3]))
        points3d = []
        for i in range(len(ir)):
            points3d.append(encodeData(pmat[:,i,:]))

        if func2:
            pmat = func2(pts).reshape(ne, len(ir)//2, 4)
            funcmin = min(funcmin, np.min(pmat))
            funcmax = max(funcmax, np.max(pmat))
            for i in range(len(ir)//2):
                points3d.append(encodeData(pmat[:,i,:]))
        d['points3d'] = points3d
        if func:
            d['show_clipping_function'] = mesh.dim==3
    if func:
        d['show_surface_function'] = True
        d['funcmin'] = funcmin
        d['funcmax'] = funcmax
    timer4.Stop()
    timer.Stop()
    return d

def Draw(mesh_or_func, mesh_or_none=None, name='function', websocket=False, *args, **kwargs): # currently assumes 2D mesh
    if 'order' in kwargs:
        order=kwargs['order']
    else:
        order = 2
    if isinstance(mesh_or_func, ngs.Mesh):
        mesh = mesh_or_func
        func = None

    if isinstance(mesh_or_func, ngs.CoefficientFunction):
        func = mesh_or_func
        mesh = mesh_or_none

    if isinstance(mesh_or_func, ngs.GridFunction):
        func = mesh_or_func
        mesh = mesh_or_none or func.space.mesh
        
    d = BuildRenderData(mesh, func, order=order)
    if websocket:
        d['websocket_url'] = "ws://localhost:" + str(websocket_port)

    scene = WebGLScene(func, mesh, order, websocket)
    writeHTML(d)
    return scene

tencode = ngs.Timer("encode")
def encodeData( array ):
    from base64 import b64encode
    tencode.Start()
    values = np.array(array.flatten(), dtype=np.float32)
    res = b64encode(values).decode("ascii")
    tencode.Stop()
    return res

