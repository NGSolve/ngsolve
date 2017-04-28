from ngsolve import *
from netgen.geom2d import SplineGeometry

geom = SplineGeometry("cavity.in2d")
# mesh = Mesh("cavity.vol.gz")
mesh = Mesh(geom.GenerateMesh (maxh=0.05))
print ("nv = ", mesh.nv)
mesh.Curve(5)
# define constant hpref = 4

SetPMLParameters (rad=0.8, alpha=2)
p=pml.Radial((0,0),0.8,2j)
mesh.SetPML(p,'pml')

fes = H1(mesh, order=4, complex=True, dirichlet=[3])
u = fes.TrialFunction()
v = fes.TestFunction()

a = BilinearForm (fes, symmetric=True)
a += BFI("PML_laplace", coef=1)

m = BilinearForm (fes, symmetric=True)
m += BFI("PML_mass", coef=1)
a.Assemble()
m.Assemble()

u = GridFunction(fes, multidim=50)
lam = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), u.vecs, shift=400)
print ("lam: ", lam)

Draw (u)

