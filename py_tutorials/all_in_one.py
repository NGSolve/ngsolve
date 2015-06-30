from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *

from numpy import linspace
from timeit import Timer

mesh = Mesh("square.vol")


v = FESpace ("h1ho", mesh, order=8, dirichlet=[1])
u = GridFunction (v)

f = LinearForm (v)
f += LFI (name = "source", coef = 1, 
          flags = { "something" : 123 }, definedon = [0,1,2])
f.Assemble()
# print (Timer (f.Assemble).timeit(number=1000))


a = BilinearForm (v, flags = { "symmetric" : True, "eliminate_internal" : True })
a += BFI ("mass", coef = 1)
a += BFI ("laplace", coef = 1)


c = Preconditioner (a, "multigrid", { "test" : True, "smoother" : "block", "finesmoothingsteps" : 1 })
# c = Preconditioner (a, "bddc", { "test" : True })

print ("now assemble the matrix ...")
a.Assemble()
c.Update()

# inv = a.mat.Inverse(v.FreeDofs())
# u.vec.data = inv * f.vec

solver = CGSolver (a.mat, c.mat, printrates=True)
u.vec.data = solver * f.vec

# choose this one if static condensation is applied:
# np = BVP(bf = a, lf = f, gf = u, pre=c)
# np.Do()

Draw (u)

sampling = [ (x,y,u(x,y)) for x in linspace(0,1,6) for y in linspace(0,1,6) ]


# just an idea by now ...
# sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) if mesh.Contains(x,y)]

