geometry = square.in2d
mesh = square.vol

# remove this line no windows:
shared = libmyngsolve


define coefficient csource
4000, 

define coefficient lam
1, 


define coefficient penalty
1e6, 1e6, 1e6, 1e6, 


define fespace v -order=3
define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v
laplace lam
mynonlinear
robin penalty

define linearform f -fespace=v
source csource

numproc nonlinearsolve np1 -bilinearform=a -linearform=f -gridfunction=u -maxit=50


