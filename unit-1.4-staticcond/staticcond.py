from netgen.geom2d import unit_square
from ngsolve import Mesh, H1, BilinearForm, grad, LinearForm
from ngsolve import SymbolicBFI, SymbolicLFI, GridFunction, Draw

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

fes = H1(mesh, order=10, dirichlet=[1, 2])
u = fes.TestFunction()
v = fes.TrialFunction()

a = BilinearForm(fes, flags={"eliminate_internal": True})
a += SymbolicBFI(grad(u) * grad(v))
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(1 * v)
f.Assemble()

u = GridFunction(fes)

f.vec.data += a.harmonic_extension_trans * f.vec
u.vec.data = a.mat.Inverse(fes.FreeDofs(True)) * f.vec
u.vec.data += a.harmonic_extension * u.vec
u.vec.data += a.inner_solve * f.vec

Draw(u)
