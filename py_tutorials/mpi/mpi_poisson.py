# mpirun -np 5 ngspy mpi_poisson.py
from ngsolve.ngstd import MPIManager
MPIManager.InitMPI()

rank = MPIManager.GetRank()
print("my rank is ", rank)

from netgen.geom2d import unit_square
if rank==0:
    mesh = unit_square.GenerateMesh(maxh=0.07)
    mesh.Save("square.vol")
MPIManager.Barrier()

# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
import netgen.meshing as netgen

# generate a triangular mesh of mesh-size 0.2
ngmesh = netgen.Mesh(dim=2)
ngmesh.Load("square.vol")
mesh = Mesh(ngmesh)

# H1-conforming finite element space
V = H1(mesh, order=3, dirichlet=[1,2,3,4])

# the right hand side
f = LinearForm (V)
f += Source (32 * (y*(1-y)+x*(1-x)))

# the bilinear-form 
a = BilinearForm (V, symmetric=False)
a += Laplace (1)

a.Assemble()
f.Assemble()

# the solution field 
u = GridFunction (V)
u.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="mumps") * f.vec

exact = 16*x*(1-x)*y*(1-y)
print ("L2-error:", sqrt (Integrate ( (u-exact)*(u-exact), mesh)))
