geometry = square.in2d
mesh = square.vol

define constant heapsize = 1000000
shared = libmyngsolve



define coefficient flow
( (y-0.5), (0.5-x) ),

define coefficient u0
(exp(-40*((x-0.7)*(x-0.7)+(y-0.7)*(y-0.7)))),
(x)

define fespace v -order=5 -type=l2ho -all_dofs_together

define gridfunction u -fespace=v 


numproc setvalues np1 -gridfunction=u -coefficient=u0


numproc linhyp np2 -gridfunction=u -flow=flow -dt=0.001 -tend=10

numproc visualization npv1 -scalarfunction=u -subdivision=3 -nolineartexture -minval=0 -maxval=1


