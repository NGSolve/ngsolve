# Call with:
# mpirun -np 5 ngspy mpi_timeDG.py

# Then make movie with:
# python3 make_timeDG_movie.py
# (this only works if you have paraview+python3 configured properly)

# circular convection; time-DG with skeleton-formulation


from netgen.geom2d import unit_square
import netgen.meshing


from ngsolve import *

SetNumThreads(1)

from mpi4py.MPI import COMM_WORLD as comm
rank = comm.rank
np = comm.size

do_vtk = False

if rank==0:
    ngmesh = unit_square.GenerateMesh(maxh=0.05)
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
mesh = Mesh(ngmesh)

# build L2-FESpace
fes = L2(mesh, order=4)

# Trial- and Test-functions
u = fes.TrialFunction()
v = fes.TestFunction()

# RHS
b = CoefficientFunction( (y-0.5,0.5-x) )
bn = b*specialcf.normal(2)

# Skeleton DG-formulation
a = BilinearForm(fes)
a += SymbolicBFI (-u * b*grad(v))
a += SymbolicBFI ( bn*IfPos(bn, u, u.Other()) * (v-v.Other()), VOL, skeleton=True)
a += SymbolicBFI ( bn*IfPos(bn, u, 0) * v, BND, skeleton=True)

u = GridFunction(fes)
u.Set(exp (-40 * ( (x-0.7)*(x-0.7) + (y-0.7)*(y-0.7) )))

w = u.vec.CreateVector()

t = 0
tau = 5e-4
tend = 2
count = 0

vtk_interval = int(0.02/tau);

if do_vtk:
    import os
    output_path = os.path.dirname(os.path.realpath(__file__)) + "/timeDG_output"
    if rank==0 and not os.path.exists(output_path):
        os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    vtk = VTKOutput(ma=mesh,coefs=[u],names=["sol"],filename=output_path+"/vtkout",subdivision=2)

with TaskManager():
    while t < tend:
        if rank==0:
            print("t = ", t)
        a.Apply (u.vec, w)
        fes.SolveM (rho=CoefficientFunction(1), vec=w)
        u.vec.data -= tau * w
        t += tau

        if count%vtk_interval==0 and do_vtk:
            vtk.Do(time=t)
        count = count+1;

        comm.Barrier()
