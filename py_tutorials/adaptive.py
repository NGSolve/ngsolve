from ngsolve import *
from netgen.occ import *


#   sub-domains "base", "top", "chip", Dirichlet boundary "bot"
#
#
#             +-------------+
#             |             |
#             |     top     |
#             |             |
#      +------+-------------+------+
#      |                           |
#      |             +             |
#      |           /   \           |
#      |          +  chip +        |
#      |           \   /    base   |
#      |             +             |
#      |                           |
#      +---------------------------+
#                   bot
#

def MakeGeometry():
    base = Rectangle(1, 0.6).Face()
    chip = MoveTo(0.5,0.15).Line(0.15,0.15).Line(-0.15,0.15).Line(-0.15,-0.15).Close().Face()
    top = MoveTo(0.2,0.6).Rectangle(0.6,0.2).Face()
    base -= chip

    base.faces.name = "base"
    chip.faces.name = "chip"
    top.faces.name = "top"

    # glue faces to one geometry with conforming interfaces
    geo = Glue([base,chip,top])
    geo.edges.name = "default"
    geo.edges.Min(Y).name = "bot"
    return geo



mesh = MakeGeometry().GenerateMesh(dim=2, maxh=0.2)


fes = H1(mesh, order=3, dirichlet="bot")
u,v = fes.TnT()

# one heat conductivity coefficient per sub-domain
lam = mesh.MaterialCF({ "base" : 1, "top" : 1000, "chip" : 10 })
a = BilinearForm(lam*grad(u)*grad(v)*dx)

# heat-source in the chip
f = LinearForm(1*v*dx("chip"))

c = preconditioners.MultiGrid(a, inverse = "sparsecholesky")

gfu = GridFunction(fes)
Draw (gfu)

# finite element space and gridfunction to represent
# the heatflux:
space_flux = HDiv(mesh, order=2)
gf_flux = GridFunction(space_flux, "flux")

def SolveBVP():
    a.Assemble()
    f.Assemble()
    inv = CGSolver(a.mat, c.mat)
    gfu.vec.data = inv * f.vec
    Redraw (blocking=True)



l = []

def CalcError():
    flux = lam * grad(gfu)
    # interpolate finite element flux into H(div) space:
    gf_flux.Set (flux)

    # Gradient-recovery error estimator
    err = 1/lam*(flux-gf_flux)*(flux-gf_flux)
    elerr = Integrate (err, mesh, VOL, element_wise=True)

    maxerr = max(elerr)
    l.append ( (fes.ndof, sqrt(sum(elerr)) ))
    print ("maxerr = ", maxerr)

    # mark elements for refinement
    mesh.ngmesh.Elements2D().NumPy()["refine"] = elerr.NumPy() > 0.25*maxerr


with TaskManager():
    while fes.ndof < 100000:  
        SolveBVP()
        CalcError()
        mesh.Refine()
    
SolveBVP()




## import matplotlib.pyplot as plt

## plt.yscale('log')
## plt.xscale('log')
## plt.xlabel("ndof")
## plt.ylabel("H1 error-estimate")
## ndof,err = zip(*l)
## plt.plot(ndof,err, "-*")

## plt.ion()
## plt.show()

## input("<press enter to quit>")

