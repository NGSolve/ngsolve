geometry = cube.geo
mesh = cube.vol

define constant heapsize = 100000000
# define constant geometryorder = 1
# define constant refinep = 1 

define constant pen = 1e7

define coefficient lam
1, 

define coefficient rho
1,

define coefficient penalty
(pen), 0,

define coefficient coef_source
(x*x*x*x), 0

define coefficient dirich_bc
(pen*x*y), 0, 



define fespace v -order=6 -type=h1ho
define fespace vp -order=5 -dim=3 -type=h1ho

define gridfunction u -fespace=v -nested


define linearform f -fespace=v
source coef_source
neumann dirich_bc

define bilinearform a -fespace=v -symmetric -linearform=f -eliminate_internal
laplace lam
mass rho
robin penalty



# define preconditioner c -type=direct -bilinearform=a
# define preconditioner c -type=local -bilinearform=a
define preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=1 -smoother=block -notest -blocktype=9
# define preconditioner c -type=amg -bilinearform=a -coefe=lam -notiming -test

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=200 -noprint -prec=1e-8
#numproc calcflux np2 -bilinearform=a -solution=u -flux=p -applyd
#numproc drawflux np4 -order
#numproc drawflux np3 -bilinearform=aid -solution=u -noapplyd -label=sol
        