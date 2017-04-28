
from netgen.csg import *

geo = CSGeometry()

box = OrthoBrick(Pnt(0,0,0),Pnt(1,1,1))
top = Plane(Pnt(0,0,0.52),Vec(0,0,1))
bot = Plane(Pnt(0,0,0.48),Vec(0,0,-1))
plate = box * top * bot

geo.Add((box-top).mat("air"))
geo.Add(plate.mat("plate"))
geo.Add((box-bot).mat("air"))

slices = [2**(-i) for i in reversed(range(1,6))]
geo.CloseSurfaces(bot,top,slices)
nmesh = geo.GenerateMesh(maxh=0.3)
ZRefinement(nmesh,geo)

from ngsolve import *
mesh = Mesh(nmesh)
Draw(mesh)
