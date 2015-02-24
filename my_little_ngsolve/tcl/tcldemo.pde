mesh = square.vol

shared = libtcldemo

define fespace v -type=h1ho -order=5
define gridfunction u -fespace=v

numproc tcldemo np1

numproc visualization npd -scalarfunction=u -subdivision=3

