geometry = chip.in2d
mesh = chip.vol

define coefficient lam
1, 1000, 10,

define coefficient penalty
1e6, 0, 

define coefficient coef_source
0, 0, 1,

define fespace v -h1ho -order=3
# define fespace v -order=1
define fespace verr -l2 -order=0

define gridfunction u -fespace=v -nested
define gridfunction err -fespace=verr

define bilinearform a -fespace=v -symmetric
laplace lam
robin penalty

define linearform f -fespace=v
source coef_source

define preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=3

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=1000

numproc drawflux np2 -bilinearform=a -solution=u -applyd -label=flux 

numproc zzerrorestimator np3 -bilinearform=a -linearform=f -solution=u -error=err -minlevel=1
numproc markelements np4 -error=err -minlevel=1 -factor=0.5 



