geometry = square.in2d
mesh = square.vol

shared = libmixed


define coefficient coefa
1,

define coefficient coef_source
1,

define coefficient one
1,


define fespace v -type=hybridmixed -order=3 -dirichlet=[1,2]

define gridfunction u -fespace=v -nested

define linearform f -fespace=v
source coef_source -comp=2

define bilinearform a -fespace=v -symmetric -eliminate_internal -keep_internal  -linearform=f
hybridmixeddiffusion coefa

define preconditioner c -type=direct -bilinearform=a -inverse=pardiso

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u  -preconditioner=c -print


define bilinearform avisu -fespace=v -nonassemble
mass one -comp=2
numproc drawflux npdf1 -solution=u -bilinearform=avisu -label=u

define bilinearform avissigma -fespace=v -nonassemble
masshdiv one -comp=1
numproc drawflux npdf2 -solution=u -bilinearform=avissigma -label=sigma



