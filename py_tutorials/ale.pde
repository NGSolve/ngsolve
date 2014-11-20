mesh = square.vol

pymodule = ale

constant numthreads = 1

fespace v -type=h1ho -order=3 -dirichlet=[1,2]
fespace vd -type=h1ho -order=3 -vec

gridfunction u -fespace=v
gridfunction def -fespace=vd

bilinearform a -fespace=v -symmetric
laplace 1

linearform f -fespace=v
source 1


pynumproc npALE np1



