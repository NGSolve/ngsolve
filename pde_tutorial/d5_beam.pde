#
# PDE example file for solving a linear elasticity model
#


# load prepared geometry and mesh file:
#
geometry = beam.geo
mesh = beam.vol


define constant heapsize = 100000000


# the elastic modulus for each sub-domain (here just one)
define coefficient E
2.1E11,

# the Poisson ratio
define coefficient nu
0.2,


# Displacement boundary conditions (clumped) are set by penalty
# Give a large value for each clamped piece of the boundary 
# The index corresponds to the value set in 'Edit Boundary Conditons'

define coefficient penalty
1e20, 0, 


# A volume force in z-direction
define coefficient coef_force_z
7e4, 

# A surface force in z-direction for each piece of the boundary
define coefficient coef_surface_force_z
0, 1e5, 


define fespace v -type=h1ho -dim=3 -order=4 -dirichlet=[1]
define fespace vp -type=h1ho -dim=6  -order=3

define gridfunction u -fespace=v
define gridfunction stress -fespace=vp

# generate load vector, volume force in z-direction (-comp=3)
define linearform f -fespace=v
# source coef_force_z -comp=1
neumann coef_surface_force_z -comp=3

# define system matrix. robin adds penalty terms to the x,y, and z-components
define bilinearform a -fespace=v -symmetric  -eliminate_internal -linearform=f
elasticity E nu


# use either a direct factorization, or a multigrid preconditioner.
# for problems smaller than 10000 nodes, the direct one is usually faster,
# otherwise the multigrid is faster (and needs much less memory)
 
# define preconditioner c -type=direct -bilinearform=a
define preconditioner c -type=multigrid -bilinearform=a  # -inverse=sparsecholesky

#solve the problem by calling the cg-iteration
numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=200


# postprocessing computes the stresses:
numproc calcflux np2 -bilinearform=a -solution=u -flux=stress -applyd




numproc visualization npv1 -vectorfunction=u -subdivision=2 -nolineartexture -deformationscale=100 -scalarfunction=stress.1
