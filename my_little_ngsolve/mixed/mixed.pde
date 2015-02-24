geometry = square.in2d
mesh = square.vol

shared = libmixed


define coefficient coefa
1,

define coefficient coef_source
1,

define coefficient one
1,

define coefficient reg
-1e-8,



define fespace v -type=mixed -order=3

define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v -symmetric 
mixeddiffusion coefa
mass reg -comp=2

define linearform f -fespace=v
source coef_source -comp=2



numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u  -solver=direct


define bilinearform avisu -fespace=v -nonassemble
mass one -comp=2
numproc drawflux npdf1 -solution=u -bilinearform=avisu -label=u

define bilinearform avissigma -fespace=v -nonassemble
masshdiv one -comp=1
numproc drawflux npdf2 -solution=u -bilinearform=avissigma -label=sigma


# numproc shapetester nptest -gridfunction=u
