from netgen.geom2d import *
from ngsolve import *

geo = SplineGeometry()
geo.AddRectangle ( (0,0), (2,2), bcs=["b","r","t","l"], leftdomain=1)
geo.AddRectangle ( (1,0.9), (1.3,1.4), bcs=["b2","r2","t2","l2"], leftdomain=2, rightdomain=1)
geo.SetMaterial (1, "outer")
geo.SetMaterial (2, "inner")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

fes1 = H1(mesh, definedon="inner")
u1 = GridFunction(fes1, "u1")
u1.Set (x*y)

fes = H1(mesh, order=3,dirichlet="b|l|r")
u = fes.TrialFunction()
v = fes.TestFunction()

gfu = GridFunction(fes)

f = LinearForm(fes)
f += SymbolicLFI (u1*v, definedon=mesh.Materials("inner"))
f.Assemble()

a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v))
a.Assemble()

gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw (u1)
Draw (gfu)
