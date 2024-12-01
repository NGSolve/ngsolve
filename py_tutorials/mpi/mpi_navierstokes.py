import netgen.meshing
from netgen.geom2d import SplineGeometry


from ngsolve import *

SetNumThreads(1)  

from mpi4py.MPI import COMM_WORLD as comm
rank = comm.rank
np = comm.size

do_vtk = False

# viscosity
nu = 0.001

# timestepping parameters
tau = 0.001
tend = 3

# mesh = Mesh("cylinder.vol")
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (2, 0.41), bcs = ("wall", "outlet", "wall", "inlet"))
geo.AddCircle ( (0.2, 0.2), r=0.05, leftdomain=0, rightdomain=1, bc="cyl")
if rank==0:
    ngmesh = geo.GenerateMesh(maxh=0.08)
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
    ngmesh.SetGeometry(geo)
    
mesh = Mesh(ngmesh)

#does not work with mpi yet...
#mesh.Curve(3)

V = H1(mesh,order=3, dirichlet="wall|cyl|inlet")
Q = H1(mesh,order=2)

X = FESpace([V,V,Q])

ux,uy,p = X.TrialFunction()
vx,vy,q = X.TestFunction()

div_u = grad(ux)[0]+grad(uy)[1]
div_v = grad(vx)[0]+grad(vy)[1]

stokes = nu*grad(ux)*grad(vx)+nu*grad(uy)*grad(vy)+div_u*q+div_v*p - 1e-10*p*q
a = BilinearForm(X)
a += SymbolicBFI(stokes)
a.Assemble()

# nothing here ...
f = LinearForm(X)   
f.Assemble()

# gridfunction for the solution
gfu = GridFunction(X)

# parabolic inflow at bc=1:
uin = 1.5*4*y*(0.41-y)/(0.41*0.41)
gfu.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))

velocity = CoefficientFunction(gfu.components[0:2])


# solve Stokes problem for initial conditions:
#inv_stokes = a.mat.Inverse(X.FreeDofs(), inverse="mumps")
inv_stokes = a.mat.Inverse(X.FreeDofs(), inverse="masterinverse")
res = f.vec.CreateVector()
res.data = f.vec - a.mat*gfu.vec
gfu.vec.data += inv_stokes * res

# matrix for implicit Euler 
mstar = BilinearForm(X)
mstar += SymbolicBFI(ux*vx+uy*vy + tau*stokes)
mstar.Assemble()

# inv = mstar.mat.Inverse(X.FreeDofs(), inverse="masterinverse")
inv = mstar.mat.Inverse(X.FreeDofs(), inverse="masterinverse")

# the non-linear term 
conv = BilinearForm(X, nonassemble = True)
conv += SymbolicBFI( CoefficientFunction( (ux,uy) ) * (grad(ux)*vx+grad(uy)*vy) )

t = 0
vtk_interval = int(0.05/tau);

import os
#output_path = os.path.dirname(os.path.realpath(__file__)) + "/navierstokes_output"
output_path = os.path.dirname(os.path.realpath(__file__)) + "/navierstokes_output"
if rank==0 and not os.path.exists(output_path):
    os.mkdir(output_path)
comm.Barrier() #wait until master has created the directory!!

if do_vtk:
    vtk = VTKOutput(ma=mesh,coefs=[velocity],names=["u"],filename=output_path+"/vtkout",subdivision=2)
    vtk.Do()

count = 1;
# implicit Euler/explicit Euler splitting method:
with TaskManager():
    while t < tend:
        if rank==0:
            print ("t=", t)
            
        conv.Apply (gfu.vec, res)
        res.data += a.mat*gfu.vec
        gfu.vec.data -= tau * inv * res

        if count%vtk_interval==0 and do_vtk:
            vtk.Do(time = t)
        count = count+1;

        t = t + tau
        comm.Barrier()
