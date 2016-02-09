#
# geometric non-linear elasticity with Neo-Hook hyperelastic material
#
# featuring automatic differentiation in SymbolicEnergy
#

import netgen.geom2d as geom2d
from ngsolve import *

geo = geom2d.SplineGeometry()
pnums = [ geo.AddPoint (x,y,maxh=0.01) for x,y in [(0,0), (1,0), (1,0.1), (0,0.1)] ]
for p1,p2,bc in [(0,1,"bot"), (1,2,"right"), (2,3,"top"), (3,0,"left")]:
     geo.Append(["line", pnums[p1], pnums[p2]], bc=bc)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))


E, nu = 210, 0.2
mu  = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

fes = H1(mesh, order=2, dirichlet="left", dim=mesh.dim)

u  = fes.TrialFunction()

force = CoefficientFunction( (0,1) )

# some utils ...
def IdentityCF(dim):
    return CoefficientFunction( tuple( [1 if i==j else 0 for i in range(dim) for j in range(dim)]), dims=(dim,dim) )

def Trace(mat):
    return sum( [mat[i,i] for i in range(mat.dims[0]) ])

def Det(mat):
    if mat.dims[0] == 2:
        return mat[0,0]*mat[1,1]-mat[0,1]*mat[1,0]

I = IdentityCF(mesh.dim)
F = I + u.Deriv()   # attention: row .. component, col .. derivative
C = F * F.trans  
E = 0.5 * (C-I)

def Pow(a, b):
    return exp (log(a)*b)
  
def NeoHook (C):
    return 0.5 * mu * (Trace(C-I) + 2*mu/lam * Pow(Det(C), -lam/2/mu) - 1)



factorspace = FESpace("number", mesh)
factor = GridFunction(factorspace)
factor.vec[0] = 0.1

  

a = BilinearForm(fes, symmetric=False)
a += SymbolicEnergy(  NeoHook (C).Compile() )
a += SymbolicEnergy(  (-factor * InnerProduct(force,u) ).Compile() )


u = GridFunction(fes)
u.vec[:] = 0

res = u.vec.CreateVector()
w = u.vec.CreateVector()


for loadstep in range(50):
    
    print ("loadstep", loadstep)
    factor.vec[:] = (loadstep+1)/10
    
    for it in range(5):
        print ("Newton iteration", it)
        print ("energy = ", a.Energy(u.vec))
        a.Apply(u.vec, res)
        a.AssembleLinearization(u.vec)
        inv = a.mat.Inverse(fes.FreeDofs() ) 
        w.data = inv*res
        print ("err^2 = ", InnerProduct (w,res))
        u.vec.data -= w
    
    Draw (u, mesh, "displacement")
    SetVisualization (deformation=True)
    input ("<press a key>")



    
