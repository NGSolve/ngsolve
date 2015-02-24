# geometry = beam.in2d
# mesh = beam.vol.gz
# geometry = plate.geo
# mesh = plate.vol.gz
 geometry = pipe.geo
 mesh = pipe.vol.gz

# define constant testout = "test.out"
# define constant numthreads = 1
define constant heapsize = 20000000
define constant geometryorder = 3
# shared = libHDG

define constant E = 1000
define constant nu = 0.2
define coefficient cf  (1,0,0)
define coefficient cfz 1

define fespace v -type=h1ho -order=3 -dirichlet=[1] -dim=3

define gridfunction u -fespace=v -addcoef

define coefficient cg 0, 0, 1
define linearform f -fespace=v
neumann cg -comp=3

define bilinearform a -fespace=v -symmetric -noprintelmat -xxelmatev
elasticity E nu 

numproc bvp np1 -gridfunction=u -bilinearform=a -linearform=f -solver=direct

