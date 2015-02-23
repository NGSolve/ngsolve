from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *

import numpy as np
from timeit import Timer

mesh = Mesh("square.vol")


v = FESpace ("h1ho", mesh, order=6, dirichlet=[1])
v.Update()

u = GridFunction (v)
u.Update()

f = LinearForm (v)
f.Add (LFI (name = "source", dim = 2, coef = ConstantCF(1),
            flags = { "something" : 123 }, definedon = [0,1,2]))
f.Assemble()
# print (Timer (f.Assemble).timeit(number=1000))

a = BilinearForm (v, flags = { "symmetric" : True, "eliminate_internal" : False })
a.Add (BFI ("mass", 2, ConstantCF(1)))
a.Add (BFI ("laplace", 2, ConstantCF(1)))


c = Preconditioner (a, "multigrid", { "test" : True, "smoother" : "block" })
# c = Preconditioner (a, "bddc", { "test" : True })

print ("now assemble the matrix ...")
a.Assemble()
c.Update()

# inv = a.mat.Inverse(v.FreeDofs())
# u.vec.data = inv * f.vec

solver = CGSolver (a.mat, c.mat, printrates=True)
u.vec.data = solver * f.vec


sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) ]

# sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) if mesh.Contains(x,y)]



