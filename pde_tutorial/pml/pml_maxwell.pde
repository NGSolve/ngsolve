geometry = arperture.geo
mesh = arperture.vol.gz

define constant geometryorder = 3

define constant pml_r = 2

define constant one = 1
define constant k = (2*pi)
define constant mkk = (-k*k)

define fespace v -type=hcurlho -order=2 -complex

define gridfunction u -fespace=v


define coefficient alpha
(-k), (0), 0, 

define coefficient gin 
(0,1,0), (0,0,0), (0,0,0),

define linearform f -fespace=v
neumannedge gin -definedon=1


define bilinearform a -fespace=v -symmetric
curlcurledge one -definedon=1
massedge mkk     -definedon=1
PML_curlcurledge one -definedon=2
PML_massedge mkk     -definedon=2
robinedge alpha -imag

define preconditioner c -type=direct -bilinearform=a # -inverse=pardiso

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c


