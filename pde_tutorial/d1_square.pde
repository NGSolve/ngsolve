#
# solve the Poisson equation -Delta u = f
#
# with boundary conditions
#      u = 0  on Gamma1 and Gamma2
#  du/dn = 0  on Gamma3 and Gamma4


# load geometry
geometry = square.in2d

# and mesh
mesh = square.vol


# coefficient for Laplae
define coefficient lam
1,

#
# Dirichlet boundary conditions are implemented by penalty
# large value on Gamma1 and Gamma2, on penalty on Gamma3 and Gamma4
#
define coefficient penalty
1e5, 0, 0, 0, 

# coefficient for the source
define coefficient coef_source
1,


# define a second order fespace (play around with -order=...)
define fespace v -order=5

# the solution field ...
define gridfunction u -fespace=v -nested

# the bilinear-form. 
define bilinearform a -fespace=v -symmetric
laplace lam
robin penalty

define linearform f -fespace=v
source coef_source

define preconditioner c -type=direct -bilinearform=a
# define preconditioner c -type=local -bilinearform=a
# define preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=1 


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000

define bilinearform ag -fespace=v -nonassemble
laplace lam

numproc drawflux np2 -bilinearform=ag -solution=u -label=flux -applyd


numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture