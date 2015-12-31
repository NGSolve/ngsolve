########################################################################################
# solve poisson equation
#    - eps Delta u = f in Omega,  u=g on Gamma0
########################################################################################


geometry = square.in2d
mesh = square.vol

# geometry = cube.geo
# mesh = cube.vol


define constant heapsize = 10000000


# ho-fespace: compound space for (u, lambda) ######
define fespace v -type=HDG -order=3 -dirichlet=[1,2]


define gridfunction u -fespace=v




## boundary terms ########################################

define coefficient coef_dirichlet
0, (y * (1-y)), 0, 0,  

numproc setvalues npsv -coefficient=coef_dirichlet -gridfunction=u -boundary -comp=2


## some coefficients #####################################
define coefficient one
1, 

define coefficient alpha
2,


define coefficient b
(50,10),

define coefficient cf 
1,


define linearform f -fespace=v
# source cf -comp=1


define bilinearform a -fespace=v -eliminate_internal -linearform=f -printelmat
HDG_laplace one alpha
HDG_convection b


define preconditioner c -type=direct -bilinearform=a # -inverse=pardiso

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c




# visualization

## automatic visualizaiton for 2D
numproc visualization vis1 -scalarfunction=u -subdivision=2

## automatic visualization for 3D
# numproc visualization vis2 -scalarfunction=u -subdivision=1 -clipvec=[0,0,-1] -clipsolution=scalar

