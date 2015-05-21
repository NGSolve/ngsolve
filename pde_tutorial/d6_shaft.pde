geometry = shaft.geo
mesh = shaft.vol

define constant geometryorder = 3

define coefficient E
1e6, 

define coefficient nu
0, 

define coefficient transmissionbc
0, 0, 0, 0, 0.01, 0, 0, 0, 0,

define coefficient bearingbc1
1, 0, 0, 0, 0, 0, 0,

define coefficient bearingbc2
0, 1, 0, 0, 0, 0, 0,

define coefficient bearingbc3
0, 0, 1, 0, 0, 0, 0, 

define coefficient surfaceloadz
0, 0, 0, 1, 0, 0, 0, 0, 0,


define fespace v -type=h1ho -dim=3 -order=3 
define fespace vp -type=h1ho -dim=6 -order=2
define fespace verr -type=l2 -order=0

define gridfunction u -fespace=v -nested #-addcoef
define gridfunction p -fespace=vp
define gridfunction error -fespace=verr

define bilinearform a -fespace=v -symmetric
elasticity E nu
robin transmissionbc -comp=1
robin transmissionbc -comp=2
robin transmissionbc -comp=3

define linearform f -fespace=v
neumann surfaceloadz -comp=3

define linearform cnsty1 -fespace=v
neumann bearingbc1 -comp=2
define linearform cnstz1 -fespace=v
neumann bearingbc1 -comp=3

define linearform cnsty2 -fespace=v
neumann bearingbc2 -comp=2
define linearform cnstz2 -fespace=v
neumann bearingbc2 -comp=3

define linearform cnsty3 -fespace=v
neumann bearingbc3 -comp=2
define linearform cnstz3 -fespace=v
neumann bearingbc3 -comp=3



define preconditioner c -type=multigrid -bilinearform=a -smoothingsteps=2
# define preconditioner c -type=local -bilinearform=a
# define preconditioner c -type=direct -bilinearform=a


numproc constrainedbvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=50 -constraints=[cnsty1,cnstz1,cnsty2,cnstz2,cnsty3,cnstz3] 

numproc calcflux np2 -bilinearform=a -solution=u -flux=p -applyd
numproc drawflux np2a -bilinearform=a -solution=u -label=stress -applyd

numproc zzerrorestimator np3 -bilinearform=a -linearform=f -solution=u -error=error -flux=p
numproc markelements np3a -error=error -factor=0.5 -minlevel=1 
