geometry = cube.geo
mesh = cube.vol

constant heapsize = 10000000

constant pen = 1e7

coefficient lam
1, 

coefficient rho
1,

coefficient penalty
(pen), 0,

coefficient dirich_bc
(pen*x*y), 0, 

fespace v -order=4 -type=h1ho
fespace vp -order=3 -dim=3 -type=h1ho

gridfunction u -fespace=v -nested
gridfunction p -fespace=vp

linearform f -fespace=v
source x*x*x*x
neumann dirich_bc

bilinearform a -fespace=v -symmetric -linearform=f -eliminate_internal
laplace lam
mass rho
robin penalty



# preconditioner c -type=direct -bilinearform=a
# preconditioner c -type=local -bilinearform=a
preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=1 -smoother=block -notest -blocktype=9
# preconditioner c -type=amg -bilinearform=a -coefe=lam -notiming -test

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=200 -noprint -prec=1e-8

numproc calcflux np2 -bilinearform=a -solution=u -flux=p -applyd
        
numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture
