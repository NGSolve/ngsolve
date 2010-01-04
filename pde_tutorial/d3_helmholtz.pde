mesh = piezo2d40round4.vol

define constant k = 1.8

define coefficient coef_lam
1, 

define coefficient coef_mass
(-k*k),

define coefficient absorb
0, 0, 0, (-k*I), 

define coefficient coef_source
1, -1, 0, 0,

define coefficient coef_dirichlet
1e5, (-1e5), 0, 0,


define fespace v -order=3 -complex -dirichlet=[1,2]
define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v
laplace coef_lam
mass coef_mass
robin absorb

define linearform f -fespace=v
neumann coef_source


define preconditioner c -type=direct -bilinearform=a

numproc setvalues npsv -coefficient=coef_dirichlet -gridfunction=u -boundary
numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -qmr -maxsteps=5000




numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture
