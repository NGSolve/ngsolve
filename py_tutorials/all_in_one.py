from ngsolve.fem import *
from ngsolve.comp import *

import numpy as np

mesh = Mesh("square.vol")


v = FESpace ("h1ho", mesh, order=6, dirichlet=[1])
v.Update()

u = GridFunction (v)
u.Update()

f = LinearForm (v)
f.Add (LFI ("source", 2, ConstantCF(1), definedon=[0,1,2] ))
f.Assemble()

a = BilinearForm (v, flags = { "symmetric" : True, "eliminate_internal" : False })
a.Add (BFI ("mass", 2, ConstantCF(1)))
a.Add (BFI ("laplace", 2, ConstantCF(1)))


c = Preconditioner (a, "multigrid", { "test" : True, "smootherxx" : "block" })
# c = Preconditioner (a, "bddc", { "test" : True, "smoother" : "block" })

a.Assemble()
c.Update()

inv = a.mat.Inverse(v.FreeDofs())
u.vec.data = inv * f.vec

sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) ]

# sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) if mesh.Contains(x,y)]



