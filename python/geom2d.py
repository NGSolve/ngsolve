from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libgeom2d.geom2d import *
#    import libcsgvis.csgvis as csgvis
#    from libcsgvis.csgvis import MouseMove
    from libmesh.meshing import *
if __platform.startswith('win'):
    # Windows
    from nglib.geom2d import *
#    import nglib.csgvis as csgvis
#    from nglib.csgvis import MouseMove
    from nglib.meshing import *



tmp_generate_mesh = SplineGeometry.GenerateMesh

def geom2d_meshing_func (geom, **args):
    if "mp" in args:
        return tmp_generate_mesh (geom, args["mp"])
    else:
        return tmp_generate_mesh (geom, MeshingParameters (**args))


SplineGeometry.GenerateMesh = geom2d_meshing_func




unit_square = SplineGeometry()
pnts = [ (0,0), (1,0), (1,1), (0,1) ]
lines = [ (0,1,1), (1,2,2), (2,3,3), (3,0,4) ]
pnums = [unit_square.AppendPoint(*p) for p in pnts]
for l1,l2,bc in lines:
    unit_square.Append( ["line", pnums[l1], pnums[l2]], bc=bc)

all = ['SplineGeometry', 'unit_square']




