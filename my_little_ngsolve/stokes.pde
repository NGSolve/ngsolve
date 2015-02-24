geometry = square.in2d
mesh = square.vol

# remove this line no windows:
shared = libmyngsolve


define coefficient nu
1,

define coefficient reg
-1e-10,

define coefficient one
1,

define coefficient penalty
1e8, 1e8,

define coefficient ubound
0, (1e8*y*(1-y)), 

define fespace v -type=stokes  -order=5

define gridfunction up -fespace=v -addcoef

define bilinearform a -fespace=v  -symmetric
stokes nu 
mass reg -comp=3
robin penalty -comp=1
robin penalty -comp=2


define linearform f -fespace=v
neumann ubound -comp=2


define preconditioner c -type=direct -bilinearform=a


numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=up -preconditioner=c -maxsteps=100 -qmr -prec=1e-8 -print



define bilinearform evalp -fespace=v -symmetric -nonassemble 
mass one -comp=3 

# numproc drawflux np2 -bilinearform=evalu -solution=up -label=velocity
numproc drawflux np2b -bilinearform=evalp -solution=up -label=pressure


# select "vector function = velocity", and
# switch on "Draw Surface Solution"

define coefficient velocity ( (up.1, up.2) )
numproc draw npd1 -coefficient=velocity -label=velocity


numproc visualization npv1 -vectorfunction=velocity -minval=0 -maxval=0.25




