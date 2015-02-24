geometry = square.in2d
mesh = square.vol

shared = libmyngsolve



define fespace v -type=myfespace -dirichlet=[1]

define gridfunction u -fespace=v 



numproc myassembling np1


