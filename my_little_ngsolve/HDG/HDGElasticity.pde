# geometry = beam.in2d
# mesh = beam.vol.gz
# geometry = plate.geo
# mesh = plate.vol.gz
geometry = pipe.geo
mesh = pipe.vol.gz

define constant heapsize = 20000000
define constant geometryorder = 4

# constant testout = "test.out"
# constant numthreads = 1 

shared = libHDG

define constant E = 1000
define constant nu = 0.2
define coefficient cf  (0,0,1)

define fespace v -type=HDGElasticity -order=4 -dirichlet=[1] -xxprint

define gridfunction u -fespace=v -addcoef

# define coefficient cg 0, 0, 1, 
# define coefficient cg 0, 0, (2*z),
# define coefficient cg 0, 0, (1),
define coefficient ct  (0,0,0), (0,0,0), (0,0,1)


define linearform f -fespace=v
# sourceedge cf -comp=1
# neumannhdiv cg -comp=2
neumannedge ct -comp=1 -domain=3

define bilinearform a -fespace=v -symmetric -xxprintelmat -xxelmatev -eliminate_internal -xkeep_internal -xxlinearform=f -spd
HDG_elasticity E nu

define preconditioner c -type=direct -bilinearform=a


numproc bvp np1 -gridfunction=u -bilinearform=a -linearform=f -preconditioner=c

numproc drawflux npdf -solution=u -bilinearform=a -label=strain

