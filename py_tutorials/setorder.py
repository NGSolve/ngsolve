from netgen.geom2d import unit_square
from ngsolve import *

mesh = Mesh (unit_square.GenerateMesh(maxh=0.3))

fes = H1(mesh, order=1)

# defines an order policy
fes.order[EDGE] = 1     # anyway 
fes.order[TRIG] = 3
fes.order[QUAD] = 2


# sets the order for particular nodes (not yet functional)
for el in fes.Elements(VOL):
    if 0 in el.vertices:
        for e in el.edges:
            fes.order[EDGE,e] = 2
# fes.UpdateDofTable()   # missing

for el in fes.Elements(VOL):
    print (el.dofs)
    
