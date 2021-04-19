from ngsolve import *
from netgen.geom2d import SplineGeometry

geo = SplineGeometry()
geo.AddRectangle( (-0.1, -0.25), (0.1, 0.25), leftdomain=0, rightdomain=1, bc = "scatterer")
geo.AddCircle ( (0, 0), r=1, leftdomain=1, rightdomain=2)
geo.AddCircle ( (0, 0), r=1.4, leftdomain=2, rightdomain=0)

ngmesh = geo.GenerateMesh(maxh=0.04)
# ngmesh.Save("scattering.vol")
mesh = Mesh(ngmesh)
# mesh = Mesh ("scattering.vol")

mesh.SetPML(pml.Radial(origin=(0,0), rad=1, alpha=0.1j), definedon=2)

kx = 50
ky = 20
k = sqrt(kx*kx+ky*ky)

uin = exp (1J*kx*x+1J*ky*y)

fes = H1(mesh, complex=True, order=5, dirichlet="scatterer")
u = fes.TrialFunction()
v = fes.TestFunction()

uscat = GridFunction (fes)
uscat.Set (uin, definedon=mesh.Boundaries("scatterer"))

a = BilinearForm (fes, symmetric=True)
a += grad(u)*grad(v)*dx
a += -k*k*u*v*dx

f = LinearForm (fes)

a.Assemble()
f.Assemble()

res = uscat.vec.CreateVector()
res.data = f.vec - a.mat * uscat.vec
uscat.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * res


Draw (uin, mesh, "uin")
Draw (uscat, mesh, "uscat")
Draw (uin-uscat, mesh, "utot")




