########################################################################################
# solve poisson equation
#    - eps Delta u = f in Omega,  u=g on Gamma0
# by means of a mixed-hybrid method (see demo_habrid.cpp for details) 
########################################################################################


geometry = square.in2d
mesh = square.vol

# geometry = cube.geo
# mesh = cube.vol


define constant heapsize = 10000000


# ho-fespace: compound space for (u, lambda) ######
define fespace v -type=HDG -order=3 -noeliminate_internal


define gridfunction u -fespace=v


## boundary terms ########################################
define coefficient crob
1e10, 1e10, 

define coefficient cneu
0, (1e10*y * (1-y)), 0, 0,  


## some coefficients #####################################
define coefficient one
1, 

define coefficient alpha
2,


define coefficient bx
500,

define coefficient by
100,




define coefficient cf 
1,


define linearform f -fespace=v
neumann cneu -comp=2
# source cf -comp=1


define bilinearform a -fespace=v -noeliminate_internal -linearform=f -printelmat
HDG_laplace one alpha
HDG_convection bx by
robin crob -comp=2


define preconditioner c -type=direct -bilinearform=a -inverse=mumps

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -preconditioner=c





# visualization

define bilinearform drawu -fespace=v -nonassemble
mass one -comp=1

# numproc drawflux df1 -bilinearform=drawsigma -solution=u -label=sigma
numproc drawflux df2 -bilinearform=drawu -solution=u -label=u

## automatic visualizaiton for 2D
numproc visualization vis1 -scalarfunction=u -subdivision=2

## automatic visualization for 3D
# numproc visualization vis2 -scalarfunction=u -subdivision=1 -clipvec=[0,0,-1] -clipsolution=scalar

