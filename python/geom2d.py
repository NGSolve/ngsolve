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


def MakeRectangle (geo, p1, p2, **args):
    p1x, p1y = p1
    p2x, p2y = p2
    p1x,p2x = min(p1x,p2x), max(p1x, p2x)
    p1y,p2y = min(p1y,p2y), max(p1y, p2y)
    pts = [geo.AppendPoint(*p) for p in [(p1x,p1y), (p2x, p1y), (p2x, p2y), (p1x, p2y)]]
    for p1,p2 in [(0,1), (1, 2), (2, 3), (3, 0)]:
        geo.Append( ["line", pts[p1], pts[p2]], **args)

def MakeCircle (geo, c, r, **args):
    cx,cy = c
    pts = [geo.AppendPoint(*p) for p in [(cx,cy-r), (cx+r,cy-r), (cx+r,cy), (cx+r,cy+r), \
                                         (cx,cy+r), (cx-r,cy+r), (cx-r,cy), (cx-r,cy-r)]]
    for p1,p2,p3 in [(0,1,2), (2,3,4), (4, 5, 6), (6, 7, 0)]:
        geo.Append( ["spline3", pts[p1], pts[p2], pts[p3]], **args)

    
SplineGeometry.AddCircle = lambda geo, c, r, **args : MakeCircle(geo, c, r, **args)
SplineGeometry.AddRectangle = lambda geo, p1, p2, **args : MakeRectangle(geo, p1, p2, **args)
SplineGeometry.AddSegment = SplineGeometry.Append
SplineGeometry.AddPoint = SplineGeometry.AppendPoint


__all__ = ['SplineGeometry', 'unit_square']




