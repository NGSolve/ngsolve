from ngsolve.ngstd import *
from ngsolve.comp import *
from ngsolve.fem import *
from ngsolve.utils import *
from ngsolve.solve import *

from math import sqrt


# geometry = scattering.in2d
mesh = Mesh ("scattering.vol.gz")

#old version
#SetPMLParameters (rad=1, alpha=0.2)  

#new version
p=pml.Radial(origin=(0,0),rad=1,alpha=0.2j)
mesh.SetPML(p,'pml')

kx = 50
ky = 10
k = sqrt(kx*kx+ky*ky)

uin = exp (1J*kx*x+1J*ky*y)

v = FESpace("h1ho", mesh, complex=True, order=3, dirichlet=[1])

uscat = GridFunction (v, "uscat")
uscat.Set (uin,BND)

Draw (uscat)
u = v.TrialFunction()
w = v.TestFunction()
a = BilinearForm (v, symmetric=True)

#old
#a += BFI("PML_laplace", coef=1)
#a += BFI("PML_mass", coef=-k*k)

a+=SymbolicBFI(grad(u)*grad(w))
a+=SymbolicBFI(-k*k*u*w)


f = LinearForm (v)

a.Assemble()
f.Assemble()

c = Preconditioner(a, type="direct")
c.Update()

bvp = BVP(bf=a,lf=f,gf=uscat,pre=c)
bvp.Do()

Draw (uin, mesh, "uin")
Draw (uscat-uin, mesh, "utot")
