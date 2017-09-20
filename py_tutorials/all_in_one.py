from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *
from ngsolve.solve import *
from ngsolve.utils import *

from numpy import linspace
from timeit import Timer

mesh = Mesh("square.vol")


v = FESpace ("h1ho", mesh, order=4, dirichlet=[1])
u = GridFunction (v)

f = LinearForm (v)
f += Source (coef=x*y)

f.Assemble()

a = BilinearForm (v, symmetric=True, eliminate_internal = False)
a += Mass (coef=1)
a += Laplace (coef=1)


# c = Preconditioner (a, type = "multigrid", test = True, smoother = "block", finesmoothingsteps = 1)
c = Preconditioner (a, type = "bddc", test = True)

print ("now assemble the matrix ...")
a.Assemble()
c.Update()

# inv = a.mat.Inverse(v.FreeDofs())
# u.vec.data = inv * f.vec

solver = CGSolver (a.mat, c.mat, printrates=True, precision=1e-10)
u.vec.data = solver * f.vec

# choose this one if static condensation is applied:
# BVP(bf=a, lf=f, gf=u, pre=c).Do()

Draw (u)

print ("f(u) = ", f(u))
print ("A(u,u) = ", a(u,u))



# evaluate solution on square ...

sampling = [ (x,y,u(x,y)) for x in linspace(0,1,6) for y in linspace(0,1,6) ]

# just an idea by now ...
# sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6) if mesh.Contains(x,y)]

# ... and along a line

xpnts = linspace(0,1,21)
vals = [ u(x,0.5) for x in xpnts ]

import pickle
outfile = open ("linevalues.dat", "wb")
pickle.dump (xpnts, outfile)
pickle.dump (vals, outfile)
outfile.close()





# reload data for post-processing

import pickle
infile = open("linevalues.dat", "rb")
pnts2 = pickle.load(infile)
vals2 = pickle.load(infile)

import matplotlib.pyplot as plt
plt.plot(pnts2, vals2, "-*")
plt.ion()
plt.show()

