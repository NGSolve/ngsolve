from ngsolve.ngstd import *
from ngsolve.comp import *
from ngsolve.fem import *
from ngsolve.utils import *
from ngsolve.solve import *

from math import sqrt


# geometry = scattering.in2d
mesh = Mesh ("scattering.vol.gz")

SetPMLParameters (rad=1, alpha=0.2)

kx = 50
ky = 10
k = sqrt(kx*kx+ky*ky)

uin = exp (1J*kx*x+1J*ky*y)

v = FESpace("h1ho", mesh, complex=True, order=3, dirichlet=[1])

uscat = GridFunction (v, "uscat")
uscat.Set (uin)

Draw (uscat)

a = BilinearForm (v, symmetric=True)
a += BFI("PML_laplace", coef=1)
a += BFI("PML_mass", coef=-k*k)


f = LinearForm (v)

a.Assemble()
f.Assemble()

c = Preconditioner(a, type="direct")
c.Update()

bvp = BVP(bf=a,lf=f,gf=uscat,pre=c)
bvp.Do()

Draw (uin, mesh, "uin")
Draw (uscat-uin, mesh, "utot")



