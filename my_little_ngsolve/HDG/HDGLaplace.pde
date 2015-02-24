geometry = square.in2d
mesh = square.vol.gz

# define constant testout = "test.out"
# define constant numthreads = 1

shared = libHDG

define constant lam = 10
define constant cf = 1

define fespace v -type=MyHDG -order=2 -dirichlet=[1]

define gridfunction u -fespace=v

define linearform f -fespace=v
source cf -comp=1

define bilinearform a -fespace=v -symmetric -noprintelmat
MyHDG_laplace lam

numproc bvp np1 -gridfunction=u -bilinearform=a -linearform=f -solver=direct
