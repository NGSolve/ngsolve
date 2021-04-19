from netgen.geom2d import SplineGeometry
from ngsolve import *

geo = SplineGeometry()
geo.AddCircle( (0,0), 1.4, leftdomain=2)
geo.AddCircle( (0,0), 1, leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "inner")
geo.SetMaterial(2, "pml")
mesh = Mesh(geo.GenerateMesh (maxh=0.1))

mesh.SetPML(pml.Radial(rad=1,alpha=1j,origin=(0,0)), "pml")

fes = H1(mesh, order=4, complex=True)
u = fes.TrialFunction()
v = fes.TestFunction()

omega = 10

a = BilinearForm(fes)
a += (grad(u)*grad(v)-omega*omega*u*v)*dx

f = LinearForm(fes)
f += exp(-20**2*((x-0.3)*(x-0.3)+y*y))*v*dx

a.Assemble()
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse() * f.vec

Draw(gfu)
