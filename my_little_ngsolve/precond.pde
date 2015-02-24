geometry = square.in2d
mesh = square.vol

shared = libmyngsolve


define coefficient lam
1,

define coefficient penalty
1e5, 0,

define coefficient coef_source
1,


# create an instance of our new FESpace

define fespace v -type=h1ho -order=2 -dirichlet=[1]
# define fespace v -type=myfespace -secondorder -dirichlet=[1]
# define fespace v -type=myhofespace -order=5


define gridfunction u -fespace=v -nested

# use our new laplace integrator (and standard robin integrator)
define bilinearform a -fespace=v -symmetric -xxprintelmat
mylaplace lam


define linearform f -fespace=v
mysource coef_source

define preconditioner c -type=mypreconditioner -bilinearform=a -test


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000 -solver=cg

numproc visualization npvis -scalarfunction=u -subdivision=1 -nolineartexture

numproc drawflux npflux -solution=u -bilinearform=a -label=flux

# numproc shapetester nptest -gridfunction=u
