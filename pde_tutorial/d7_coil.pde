geometry = coil.geo
mesh = coil.vol


define constant geometryorder = 4
define constant heapsize = 50000000

# nu = 1/mu with mu permeability
define coefficient nu
1, 1, 

define coefficient sigma
1e-6, 1e-6,

define coefficient cs1
( x2), 0,
define coefficient cs2
( -x1 ), 0,
define coefficient cs3
0, 0, 0

define coefficient penalty
0, 0, 1e8, 0, 0, 0




define fespace v -hcurlho -order=4 -eliminate_internal -nograds 


define gridfunction u -fespace=v  -novisual

define linearform f -fespace=v
sourceedge cs1 cs2 cs3 -definedon=1


define bilinearform a -fespace=v -symmetric -linearform=f -eliminate_internal 
curlcurledge nu 
massedge sigma -order=2
robinedge penalty


define bilinearform acurl -fespace=v -symmetric -nonassemble
curlcurledge nu 


#define preconditioner c -type=multigrid -bilinearform=a  -smoother=block
define preconditioner c -type=multigrid -bilinearform=a -cylce=1 -smoother=block -blocktype=1 -coarsetype=direct -coarsesmoothingsteps=5  -notest
#define preconditioner c -type=local -bilinearform=a -test 


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c  -prec=1.e-9

numproc drawflux np5 -bilinearform=acurl -solution=u  -label=flux


# p=4, 36 its, err = 1.65647e-11, 18.92 sec
