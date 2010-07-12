#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1
#  du/dn = 1  on Gamma2


# load geometry
geometry = square.in2d

# and mesh
mesh = square.vol


# coefficient for Laplae
define coefficient lam
1,


# coefficient for the source
define coefficient coef_source
1,

# coefficient for the Neumann b.c.
define coefficient coef_neumann
0, 1,

define coefficient coef_penalty
1e5, 1e5, 1e5, 1e5,



# define a finite element space
# Dirichlet boundary is Gamma_1 
# play around with -order=...
define fespace v -order=3 -dirichlet=[1,2]

# the solution field
define gridfunction u -fespace=v -nested

define linearform f -fespace=v
source coef_source
# neumann coef_neumann

# the bilinear-form 
define bilinearform a -fespace=v -eliminate_internal -keep_internal -symmetric -linearform=f
laplace lam
# robin coef_penalty



# define preconditioner c -type=direct -bilinearform=a
# define preconditioner c -type=local -bilinearform=a 
define preconditioner c -type=multigrid -bilinearform=a -smoother=block


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000

# for the visualization of the flux
define bilinearform ag -fespace=v -nonassemble
laplace lam

numproc drawflux np2 -bilinearform=ag -solution=u -label=flux -applyd


numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture