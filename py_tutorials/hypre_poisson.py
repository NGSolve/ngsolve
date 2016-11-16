# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.ngstd import MPIManager
import netgen.meshing as netgen
MPIManager.InitMPI()

ngsglobals.msg_level = 1

ngmesh = netgen.Mesh(dim=2)
ngmesh.Load("square.vol")
mesh = Mesh (ngmesh)

# H1-conforming finite element space
V = H1(mesh, order=3, dirichlet=[1,2,3,4])

u = V.TrialFunction()
v = V.TestFunction()

# the right hand side
f = LinearForm (V)
f += SymbolicLFI (32 * (y*(1-y)+x*(1-x)) * v)

# the bilinear-form 
a = BilinearForm (V, symmetric=False)
a += SymbolicBFI(grad(u)*grad(v))


# -------
# preconditioner "hypre" does one V-cycle of boomer-amg with one sweep of sym GS smoothing
# and up to 20 levels. See ngsolve/comp/hypre_precond.cpp
# More customization is not currently supported. 
# ------

#c = Preconditioner(type="hypre", bf=a) #Set order of FESpace to 1 if using this
c = Preconditioner(type="bddc", bf=a, flags={"usehypre":True})

a.Assemble()
f.Assemble()

# the solution field 
u = GridFunction (V)

cg = CGSolver(a.mat, c.mat, printrates=True, precision=1e-12)#, maxsteps=50)

u.vec.data = cg * f.vec
#print (u.vec)


# plot the solution (netgen-gui only)
Draw (u)
Draw (-u.Deriv(), mesh, "Flux")

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (u-exact)*(u-exact), mesh)))

