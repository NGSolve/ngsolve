geometry = coil.geo
mesh = coil.vol


define constant geometryorder = 4
# define constant heapsize = 50000000

define constant mu0 = 1.257e-7

# nu = 1/mu with mu permeability
define coefficient nu
(1/mu0), (1/mu0), 

# artificial conductivity for regularization
define coefficient sigma 1e-8*nu

# 1000 Ampere-turns
define coefficient r (sqrt(x*x+y*y))
define coefficient cs 
( 1000/16 * (y/r,-x/r,0))
(0, 0, 0),

define fespace v -type=hcurlho -order=4 -nograds -dirichlet=[1]

define gridfunction u -fespace=v  -novisual

define linearform f -fespace=v
sourceedge cs  -definedon=1


define bilinearform a -fespace=v -spd -nonsym_storage -linearform=f -eliminate_internal 
curlcurledge nu 
massedge sigma


define bilinearform acurl -fespace=v -symmetric -nonassemble
curlcurledge nu 


# define preconditioner c -type=multigrid -bilinearform=a  -smoother=block
define preconditioner c -type=bddc -bilinearform=a


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c  -prec=1.e-9

numproc drawflux np5 -bilinearform=acurl -solution=u  -label=flux



numproc visualization npv1 -scalarfunction=flux:0 -clipsolution=scalar -subdivision=2 -clipvec=[0,1,0] -nolineartexture
