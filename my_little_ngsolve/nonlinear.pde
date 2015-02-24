geometry = square.in2d
mesh = square.vol

shared = libmyngsolve

define coefficient csource
4000, 

define coefficient lam
1, 


define coefficient penalty
1e6, 1e6, 1e6, 1e6, 


define fespace v -order=3 -type=h1ho
define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v -symmetric
laplace lam
mynonlinear
robin penalty

define linearform f -fespace=v
source csource

numproc nonlinearsolve np1 -bilinearform=a -linearform=f -gridfunction=u -maxit=50

numproc visualization npvis -scalarfunction=u -subdivision=2 -nolineartexture

