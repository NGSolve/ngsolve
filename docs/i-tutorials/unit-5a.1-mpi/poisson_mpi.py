from ngsolve import *
from mpi4py import MPI
from netgen.occ import unit_square

comm = MPI.COMM_WORLD

ngmesh = unit_square.GenerateMesh(maxh=0.1, comm=comm)
    
for l in range(3):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
    
fes = H1(mesh, order=3, dirichlet=".*")
u,v = fes.TnT()

a = BilinearForm(grad(u)*grad(v)*dx)
pre = preconditioners.Local(a)
a.Assemble()

f = LinearForm(1*v*dx).Assemble()
gfu = GridFunction(fes)

inv = CGSolver(a.mat, pre.mat)
gfu.vec.data = inv*f.vec

ip = InnerProduct(gfu.vec, f.vec)
printonce("(u,f) =", ip)
# (u,f) = 0.03514425357822445

    
import pickle
netgen.meshing.SetParallelPickling(True)
pickle.dump (gfu, open("solution.pickle"+str(comm.rank), "wb"))



    
