geometry = square.in2d
mesh = square.vol

shared = libmyngsolve


define coefficient lam
1,

define coefficient rho
10,


define coefficient penalty
1e5, 0, 

define coefficient coef_source
1,

define fespace v -type=h1ho -order=5

define gridfunction u -fespace=v 

define bilinearform a -fespace=v -symmetric
laplace lam
robin penalty

define bilinearform m -fespace=v -symmetric
mass rho

define linearform f -fespace=v
source coef_source


numproc parabolic np1 -bilinearforma=a -bilinearformm=m -linearform=f -gridfunction=u -dt=0.0001 -tend=100



numproc tclmenu men0 -newmenu -menuname=para -text=Parabolic
numproc tclmenu men1 -menuname=para -text=Parabolic -fieldname=u  -minval=-0.05 -maxval=0.05


numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture