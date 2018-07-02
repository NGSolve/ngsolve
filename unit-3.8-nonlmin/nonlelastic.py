import netgen.geom2d as geom2d
from ngsolve import *
from solvenlmin import *

geo = geom2d.SplineGeometry()
pnums = [ geo.AddPoint (x,y,maxh=0.01) for x,y in [(0,0), (1,0), (1,0.1), (0,0.1)] ]
for p1,p2,bc in [(0,1,"bot"), (1,2,"right"), (2,3,"top"), (3,0,"left")]:
     geo.Append(["line", pnums[p1], pnums[p2]], bc=bc)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))

# E module and poisson number:
E, nu = 210, 0.2
# Lamé constants:
mu  = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

V = H1(mesh, order=2, dirichlet="left", dim=mesh.dim)
u  = V.TrialFunction()

#gravity:
force = CoefficientFunction( (0,-1) )

# some utils:
def IdentityCF(dim):
    return CoefficientFunction( tuple( [1 if i==j else 0 for i in range(dim) for j in range(dim)]), dims=(dim,dim) )

def Trace(mat):
    return sum( [mat[i,i] for i in range(mat.dims[0]) ])

def Det(mat):
    if mat.dims[0] == 2:
        return mat[0,0]*mat[1,1]-mat[0,1]*mat[1,0]

def Pow(a, b):
    return exp (log(a)*b)
    
def NeoHook (C):
    return 0.5 * mu * (Trace(C-I) + 2*mu/lam * Pow(Det(C), -lam/2/mu) - 1)

I = IdentityCF(mesh.dim)
F = I + u.Deriv()   # attention: row .. component, col .. derivative
C = F * F.trans  

factor = Parameter(1.0)

a = BilinearForm(V, symmetric=False)
a += SymbolicEnergy(  NeoHook (C).Compile() )
a += SymbolicEnergy(  (-factor * InnerProduct(force,u) ).Compile() )

gfu = GridFunction(V)
gfu.vec[:] = 0

Draw (gfu, mesh, "u")
SetVisualization (deformation=True)

res = gfu.vec.CreateVector()
du = gfu.vec.CreateVector()

for loadstep in range(50):
    print ("loadstep", loadstep)
    factor.Set ((loadstep+1)/10)
    SolveNonlinearMinProblem(a,gfu)
    Redraw()
