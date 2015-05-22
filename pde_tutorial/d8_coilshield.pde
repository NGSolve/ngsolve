geometry = coilshield.geo
mesh = coilshield.vol

define constant geometryorder = 4
define constant refinep = 1 
define constant heapsize = 200000000

define coefficient nu material
iron 1e-4
default 1

define coefficient kappa material
default 1e-6

define coefficient cs1
0, ( y ), 0, 0, 0, 
define coefficient cs2
0, ( -x ), 0, 0, 0, 
define coefficient cs3
0, 0, 0, 0, 0, 

define coefficient penalty
0,0,0,0,0,0,0,0,
#1e+6, 0, 0, 0, 0, 0


define fespace v -type=hcurlho -order=5  -nograds

define gridfunction u -fespace=v

define linearform f -fespace=v -noprintelvec
sourceedge cs1 cs2 cs3  -definedon=2


define bilinearform a -fespace=v -symmetric  -eliminate_internal  -linearform=f
curlcurledge nu 
massedge kappa  -order=2
robinedge penalty
	
define bilinearform acurl -fespace=v -symmetric  -nonassemble
curlcurledge nu


define preconditioner c -type=multigrid -bilinearform=a -cylce=1 -smoother=block -coarsetype=direct -coarsesmoothingsteps=5 -notest


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c -maxsteps=400 -prec=1.e-8

numproc drawflux np3 -bilinearform=acurl -solution=u -label=flux
