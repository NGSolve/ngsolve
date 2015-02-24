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

# remove this line no windows:
shared = libmyngsolve


# coefficient for Laplae
define coefficient lam
1,

#
# Dirichlet boundary conditions are implemented by penalty
# large value on Gamma1 and Gamma2, on penalty on Gamma3 and Gamma4
#
define coefficient penalty
1e5, 1e5, 0, 0, 

# coefficient for the source
define coefficient coef_source
1,

define fespace v1 -type=h1ho -order=3
define gridfunction u1 -fespace=v1

define bilinearform a1 -fespace=v1 -symmetric
laplace lam
robin penalty

define linearform f1 -fespace=v1
source coef_source

define preconditioner c1 -type=direct -bilinearform=a1



define fespace v2 -type=h1ho -order=6
define gridfunction u2 -fespace=v2

define bilinearform a2 -fespace=v2 -symmetric
laplace lam
robin penalty

define linearform f2 -fespace=v2

define preconditioner c2 -type=direct -bilinearform=a2




numproc bvp np1 -bilinearform=a1 -linearform=f1 -gridfunction=u1 -preconditioner=c1

numproc democoupling np2  -gridfunction=u1 -linearform=f2
# numproc democoupling_adv np2  -gridfunction=u1 -linearform=f2

numproc bvp np3 -bilinearform=a2 -linearform=f2 -gridfunction=u2 -preconditioner=c2

