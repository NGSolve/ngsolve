from libngpy.geom2d import *
from libngpy.meshing import *

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

__all__ = ['SplineGeometry', 'unit_square']




