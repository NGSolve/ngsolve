from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import unit_square

comm = MPI.COMM_WORLD

if comm.rank == 0:
    ngmesh = unit_square.GenerateMesh(maxh=0.02)
    ngmesh.Save("square.vol")
    mesh = Mesh(ngmesh.Distribute(comm))
else:
    mesh = Mesh(netgen.meshing.Mesh.Receive(comm))


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

# if comm.size == 1:
#     Draw (gfu)
gfu.Save("solution.sol", parallel=True)





    
