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


# define a finite element space
# Dirichlet boundary is Gamma_1 
# play around with -order=...
fespace v -type=h1ho -order=3 -dirichlet=[1]

# the solution field
gridfunction u -fespace=v -nested

linearform f -fespace=v
source x*sin(pi*y)
# neumann 1 --definedon=[2]    # Neumann on Gamma_2

# the bilinear-form 
bilinearform a -fespace=v -symmetric -eliminate_internal
laplace 1


# preconditioner c -type=direct -bilinearform=a
# preconditioner c -type=local -bilinearform=a 
# preconditioner c -type=multigrid -bilinearform=a -smoother=block
preconditioner c -type=bddc -bilinearform=a


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000

numproc drawflux np2 -bilinearform=a -solution=u -label=flux -applyd

numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture
