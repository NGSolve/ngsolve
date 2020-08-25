from mpi4py import MPI
from ngsolve import *
import numpy as np
from netgen.csg import unit_cube
from ngsolve.krylovspace import CGSolver
import petsc4py.PETSc as psc

def masterprint (*args, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print (*args)


comm = MPI.COMM_WORLD

if comm.rank == 0:
    ngmesh = unit_cube.GenerateMesh(maxh=0.1).Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)


for l in range(0):
    ngmesh.Refine()
mesh = Mesh(ngmesh)


fes = H1(mesh, order=1)
u,v = fes.TnT()
a = BilinearForm(grad(u)*grad(v)*dx+u*v*ds).Assemble()
f = LinearForm(x*v*dx).Assemble()
gfu = GridFunction(fes)
inv = CGSolver(a.mat, freedofs=fes.FreeDofs(), printing=False, maxsteps=400, tol=1e-8)
gfu.vec.data = inv * f.vec
masterprint ("ngs-dot =", InnerProduct(gfu.vec, f.vec))

pardofs = fes.ParallelDofs()
globnums, nglob = pardofs.EnumerateGlobally() 

locmat = a.mat.local_mat
val,col,ind = locmat.CSR()
ind = np.array(ind, dtype='int32')

apsc_loc = psc.Mat().createAIJ(size=(locmat.height, locmat.width), csr=(ind,col,val), comm=MPI.COMM_SELF)

IS = psc.IS().createBlock (bsize=1, indices=globnums, comm=comm)
lgmap = psc.LGMap().create (bsize=1, indices=globnums, comm=comm)

mat = psc.Mat().createPython(size=nglob, comm=comm)
mat.setType(psc.Mat.Type.IS)
mat.setLGMap(lgmap)
mat.setISLocalMat(apsc_loc)
mat.assemble()

f.vec.Cumulate()

v1, v2 = mat.createVecs()

v2loc = v2.getSubVector(IS)
v2loc.getArray()[:] = f.vec.FV()
v2.restoreSubVector(IS, v2loc)

ksp = psc.KSP()
ksp.create()
ksp.setOperators(mat)
ksp.setType(psc.KSP.Type.CG)
ksp.setNormType(psc.KSP.NormType.NORM_NATURAL)
ksp.setTolerances(rtol=1e-6, atol=0, divtol=1e16, max_it=400)
ksp.solve(v2,v1)   

masterprint ("petsc-dot = ", v1.dot(v2))
masterprint ("petsc-its =", ksp.its)


