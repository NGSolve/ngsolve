mesh = piezo2d40round4.vol

define constant k = 1.8

define coefficient coef_lam
1, 
define coefficient coef_mass
(-k*k),

define coefficient absorb
0, 0, 0, (-k),

define coefficient absorb2
0, 0, 0, (-0.5/k),

define coefficient coef_source
1e5, (-1e5), 0, 0,

define coefficient penalty
1e5, 1e5, 0, 0,

define fespace v -order=2 -complex
define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v
laplace coef_lam
mass coef_mass
robin penalty
robin absorb -imag

define linearform f -fespace=v
neumann coef_source


define preconditioner c -type=direct -bilinearform=a

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -qmr -maxsteps=5000



