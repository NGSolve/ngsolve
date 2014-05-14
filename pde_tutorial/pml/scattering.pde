geometry = scattering.in2d
mesh = scattering.vol.gz

define constant kx = 50
define constant ky = 10
define constant k = (sqrt (kx*kx+ky*ky))
define coefficient uin (exp (I*kx*x+I*ky*y))


define fespace v -type=h1ho -order=3 -complex -dirichlet=[1]

define gridfunction uscat -fespace=v -addcoef
numproc setvalues np1 -gridfunction=uscat -coefficient=uin


define bilinearform a -fespace=v -symmetric
laplace 1  --definedon=[1]
mass -k*k  --definedon=[1]
PML_laplace 1 --definedon=[2]
PML_mass -k*k --definedon=[2]

define linearform f -fespace=v


define preconditioner c -bilinearform=a -type=direct

numproc bvp np2 -bilinearform=a -linearform=f -gridfunction=uscat -preconditioner=c

define coefficient utot  uscat-uin

numproc draw npd1 -coefficient=utot -label=utot
numproc draw npd1 -coefficient=uin -label=uin


numproc visualization npv1 -scalarfunction=uscat -subdivision=1 -nolineartexture
