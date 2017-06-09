# Call with: 
# mpirun -np 5 ngspy mpi_poisson.py

# Solves -laplace(u)=f on [0,1]^3

# netgen-imports
#from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import netgen.meshing as netgen

# ngsolve-imports
#  - MPIManager provides some basic MPI functionality
#  - DISTRIBUTED/CUMULATED describes the parallel status of a
#    parallel vector
from ngsolve import *
from ngsolve.ngstd import MPIManager
from ngsolve.la import DISTRIBUTED
from ngsolve.la import CUMULATED
from ngsolve.ngstd import GlobalSum

rank = MPIManager.GetRank()
np = MPIManager.GetNP()

print("Hello from rank "+str(rank)+" of "+str(np))

if rank==0:
    # master-proc generates mesh
    mesh = unit_cube.GenerateMesh(maxh=0.1)
    # and saves it to file
    mesh.Save("some_mesh.vol")

# wait for master to be done meshing
MPIManager.Barrier()

# now load mesh from file
ngmesh = netgen.Mesh(dim=3)
ngmesh.Load("some_mesh.vol")

#refine once?
# ngmesh.Refine()

mesh = Mesh(ngmesh)

# build H1-FESpace as usual
V = H1(mesh, order=3, dirichlet=[1,2,3,4])
u = fes.TrialFunction()
v = fes.TestFunction()

print("rank "+str(rank)+" has "+str(V.ndof)+" of "+str(V.ndofglobal)+" dofs!")

# RHS does not change either!
f = LinearForm (V)
f += SymbolicLFI(32 * (y*(1-y)+x*(1-x)) * v)
f.Assemble()

# neither does the BLF!
a = BilinearForm (V, symmetric=False)
a += SymbolicBFI(grad(u)*grad(v))

# Some possible preconditioners: 
#c = Preconditioner(a, type="direct", flags={"inverse":"mumps"}) # direct solve with mumps
#c = Preconditioner(a, type="bddc", flags={"inverse":"mumps"})   # BBDC + mumps for coarse inverse 
#c = Preconditioner(a, type="hypre")                             # BoomerAMG (use only for order 1)
c = Preconditioner(a, type="bddc", flags={"usehypre":True})     # BDDC + BoomerAMG for coarse matrix

a.Assemble()

# solve the equation
u = GridFunction (V)
u.vec.data = CGSolver(a.mat, c.mat) * f.vec  # use CG-solver with preconditioner c
# u.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="mumps") * f.vec  # use MUMPS parallel inverse
# u.vec.data = a.mat.Inverse(V.FreeDofs(), inverse="masterinverse") * f.vec  # use masterinverse (master proc does all the work!)


#exact solution
exact = 16*x*(1-x)*y*(1-y)
# Integrate integrates locally!
loc_error = Integrate ( (u-exact)*(u-exact), mesh)
print("local error on proc "+str(rank)+" is: "+str(loc_error))
# Sum up local error contributions
glob_error = GlobalSum(loc_error)
if rank==0:
    print ("L2-error:", sqrt(glob_error) )


# do VTK-output
import os
output_path = os.path.dirname(os.path.realpath(__file__)) + "/poisson_output"
if rank==0 and not os.path.exists(output_path):
    os.mkdir(output_path)
MPIManager.Barrier() #wait until master has created the directory!!

vtk = VTKOutput(ma=mesh, coefs=[u], names=["sol"], filename=output_path+"/vtkout_p"+str(rank), subdivision=2)
vtk.Do()
