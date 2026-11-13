#
# geometric non-linear elasticity with Neo-Hooke hyperelastic material
#
# featuring automatic differentiation of the energy functional
#

from ngsolve import *
from netgen.occ import *

# a beam with named boundaries, the mesh is refined at the corners
shape = Rectangle(1,0.1).Face()
shape.edges.Min(X).name = "left"
shape.edges.Max(X).name = "right"
shape.edges.Min(Y).name = "bot"
shape.edges.Max(Y).name = "top"
shape.vertices.maxh = 0.01
mesh = shape.GenerateMesh(dim=2, maxh=0.05)


E, nu = 210, 0.2
mu  = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

fes = VectorH1(mesh, order=2, dirichlet="left")

u  = fes.TrialFunction()

force = CF( (0,1/2) )

I = Id(mesh.dim)
F = I + Grad(u)
C = F.trans * F

def NeoHooke (C):
    return 0.5 * mu * (Trace(C-I) + 2*mu/lam * Det(C)**(-lam/2/mu) - 1)



factor = Parameter(1)

a = BilinearForm(fes)
a += Variation (NeoHooke(C).Compile()*dx)
a += Variation ((-factor * InnerProduct(force,u) ).Compile()*dx)


gfu = GridFunction(fes)

res = gfu.vec.CreateVector()
w = gfu.vec.CreateVector()


for loadstep in range(10):
    
    print ("loadstep", loadstep)
    factor.Set (loadstep+1)
    
    for it in range(5):
        print ("Newton iteration", it)
        print ("energy = ", a.Energy(gfu.vec))
        a.Apply(gfu.vec, res)
        a.AssembleLinearization(gfu.vec)
        inv = a.mat.Inverse(fes.FreeDofs() ) 
        w.data = inv*res
        print ("err^2 = ", InnerProduct (w,res))
        gfu.vec.data -= w
    
    Draw (gfu, mesh, "displacement")
    SetVisualization (deformation=True)
    input ("<press a key>")
