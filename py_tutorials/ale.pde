# mesh = square.vol
mesh = squaref.vol.gz

pymodule = ale


fespace v -type=h1ho -order=3 -dirichlet=[1,2]
fespace vd -type=h1ho -order=3 -vec

variable t = 0

gridfunction u -fespace=v
gridfunction def -fespace=vd -addcoef
gridfunction defp -fespace=vd -addcoef

bilinearform a -fespace=v -symmetric
laplace 0.1

bilinearform m -fespace=v -symmetric
mass 1

bilinearform b -fespace=v 
convection defp -transpose

linearform f -fespace=v
source 1


numproc visualization npvis -scalarfunction=u -vectorfunction=def -deformationscale=1 -subdivision=2 -nolineartexture
#  -minval=0 -maxval=0.08

# stationary equation
# pynumproc npALE np1

# parabolic equation
pynumproc npALE_instat np1

