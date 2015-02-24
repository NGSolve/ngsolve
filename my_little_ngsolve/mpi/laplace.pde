mesh = squaref.vol.gz

shared = mydd

define coefficient lam
1,

define coefficient penalty
1e5, 0,

define coefficient coef_source
1,


	
fespace v -type=h1ho -order=1 -dirichlet=[1]

gridfunction u -fespace=v -nested

bilinearform a -fespace=v -symmetric -fespace=v
laplace 1
mass 1


linearform f -fespace=v
source x

preconditioner c -type=local -bilinearform=a -test
preconditioner cdd -type=mydd -bilinearform=a  -test


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=cdd -maxsteps=1000 -prec=1e-8

numproc visualization npvis -scalarfunction=u -subdivision=0 -nolineartexture

