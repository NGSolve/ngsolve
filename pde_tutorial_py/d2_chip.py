from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *
from math import pi

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters


geom = SplineGeometry("chip.in2d")
mp = MeshingParameters (maxh=0.1)
mesh = Mesh(geom.GenerateMesh (mp))


# one coefficient per sub-domain
lam = DomainConstantCF([1, 1000, 10])

# source in sub-domain 3
source = DomainConstantCF([0, 0, 1])

v = H1(mesh, order=2, dirichlet=[1])

u = GridFunction(v)

a = BilinearForm(v, symmetric=True)
a += Laplace(lam)

f = LinearForm(v)
f += Source(source)

# c = Preconditioner(a, type="multigrid")
c = Preconditioner(a, type="local")

bvp = BVP(bf=a, lf=f, gf=u, pre=c, maxsteps=1000)

space_flux = FESpace("hdivho", mesh, order=2)
gf_flux = GridFunction(space_flux, "flux")


def SolveBVP():
    v.Update()
    u.Update()
    a.Assemble()
    f.Assemble()
    bvp.Do()
    Draw (u)



l = []

def CalcError():
    space_flux.Update()
    gf_flux.Update()

    flux = a.Flux(u)
    gf_flux.Set (flux)

    err = (flux-gf_flux)*(flux-gf_flux)
    elerr = Integrate (err, mesh, VOL, element_wise=True)

    maxerr = max(elerr)
    l.append ( (v.ndof,sum(elerr)) )
    print ("maxerr = ", maxerr)

    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, elerr[el.nr] > 0.25*maxerr)



for level in range(20):
    SolveBVP()
    CalcError()
    mesh.Refine()
    
SolveBVP()





import matplotlib.pyplot as plt

plt.yscale('log')
plt.xscale('log')
ndof,err = zip(*l)
plt.plot(ndof,err, "-*")

plt.ion()
plt.show()

