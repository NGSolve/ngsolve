geometry = square.in2d
mesh = square.vol

pymodule = instat


define fespace v -type=h1ho -order=5 -dirichlet=[1]

define gridfunction u -fespace=v 

define bilinearform a -fespace=v -symmetric
laplace 1

define bilinearform m -fespace=v -symmetric
mass 10

define linearform f -fespace=v
source 1


numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture

pynumproc npParabolic np1



