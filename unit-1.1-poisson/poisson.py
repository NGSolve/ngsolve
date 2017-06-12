# 1. Import what you need from NGSolve and Netgen Python modules

from netgen.geom2d import unit_square
from ngsolve import Mesh, H1, BilinearForm, GridFunction, SymbolicBFI
from ngsolve import SymbolicLFI, Draw, grad, x, LinearForm

# 2. Generate an unstructured mesh

mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

# 3. Declare a finite element space:

fes = H1(mesh, order=2, dirichlet="bottom|right")

# 4. Declare test function, trial function, and grid function

u = fes.TrialFunction()  # symbolic object
v = fes.TestFunction()   # symbolic object
gfu = GridFunction(fes)  # solution

# 5. Define and assemble linear and bilinear forms:

a = BilinearForm(fes, symmetric=True)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(x*v)
f.Assemble()

# 6. Solve the system

gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
Draw(gfu)
