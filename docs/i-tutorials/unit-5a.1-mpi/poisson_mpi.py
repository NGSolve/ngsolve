from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import unit_square

comm = MPI.COMM_WORLD

if comm.rank == 0:
    ngmesh = unit_square.GenerateMesh(maxh=0.1).Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)

for l in range(3):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
    
fes = H1(mesh, order=3, dirichlet=".*")
u,v = fes.TnT()

a = BilinearForm(grad(u)*grad(v)*dx)
pre = Preconditioner(a, "local")
a.Assemble()

f = LinearForm(1*v*dx).Assemble()
gfu = GridFunction(fes)


inv = CGSolver(a.mat, pre.mat)
gfu.vec.data = inv*f.vec

ip = InnerProduct(gfu.vec, f.vec)
if comm.rank == 0:
    print ("(u,f) =", ip)
# (u,f) = 0.03514425357822445

    
import pickle
netgen.meshing.SetParallelPickling(True)
pickle.dump (gfu, open("solution.pickle"+str(comm.rank), "wb"))



    
