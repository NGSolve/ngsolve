geometry = doubleglazing.in2d
mesh = doubleglazing.vol

####################### DISCONTINUOUS GALERKIN DEMO #############

#### the double glazing problem: 
#### taken from Howard Elman, David Silvester, Andy Wathen: 
#### Finite Elements and Fast Iterative Solvers 
#### with applications in incompressible fluid dynamics
#### 3.1.4 Example

define coefficient lam 
(0.005),(0.005),

#for inflow
define coefficient dirichletcoef
0, 1, 

define coefficient calpha 
10,10,

#convection velocity
define coefficient b1 
(2*y*(1-x*x)),(2*y*(1-x*x)),

#convection velocity
define coefficient b2
(-2*x*(1-y*y)),(-2*x*(1-y*y)),


define fespace vdisc -type=l2ho -order=4 -dgjumps 
define gridfunction udisc -fespace=vdisc

define bilinearform adisc -fespace=vdisc
DGIP_innfac_laplace lam calpha
DGIP_bndfac_laplace lam calpha -definedon=[1,2]
DG_innfac_convection b1 b2
DG_bndfac_convection b1 b2
convection b1 b2
laplace lam

define linearform fdisc -fespace=vdisc
DGIP_bndfac_dir lam dirichletcoef calpha -definedon=2
DG_bndfac_convdir b1 b2 dirichletcoef -definedon=2

define preconditioner c -type=direct -bilinearform=adisc

numproc bvp npdisc -bilinearform=adisc -linearform=fdisc -gridfunction=udisc -preconditioner=c

numproc visualization npvis -scalarfunction=udisc -subdivision=4 -deformationscale=1
