geometry = square.in2d
mesh = square.vol

shared = libprecond

define constant testout = "test.out"

define constant lam = 1
define coefficient cf ( sin(2*pi*x) )
define fespace v -type=myh1ho -order=5 -dirichlet=[1]


define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v -symmetric -eliminate_internal
laplace lam

define linearform f -fespace=v
source cf

# a block-Jacobi preconditioner
define preconditioner c -type=local -bilinearform=a -coarsetype=direct -block -ebe -test

# a block-Gauss-Seidel with direct solve for Direct-Solver-Clusters
# define preconditioner c -type=multigrid -bilinearform=a -smoother=block -cycle=0 -coarsetype=smoothing -ebe -test

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=100 -solver=cg

