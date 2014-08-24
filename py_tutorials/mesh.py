import sys
sys.path.append("/opt/netgen/lib")

from libngspy import *
from Ngsolve import *
from Ngcomp import *
from Ngstd import *




pde = PDE()
pde.Load ("../pde_tutorial/d1_square.pde")

mesh = pde.Mesh()

for i in mesh.Elements(VOL):
    print (i)
    print (mesh[i].vertices)


