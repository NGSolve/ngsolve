mesh = piezo2d40round4.vol.gz

define constant k = 1.8

define fespace v -type=h1ho -order=3 -complex -dirichlet=[1,2]

define gridfunction u -fespace=v -nested

define bilinearform a -fespace=v -symmetric
laplace 1
mass (-k*k)
robin (-k*I) --definedon=4

define linearform f -fespace=v

define coefficient coef_dirichlet
1, -1, 0, 0,
numproc setvalues npsv -coefficient=coef_dirichlet -gridfunction=u -boundary

numproc bvp np1 -bilinearform=a -linearform=f -gridfunction=u -solver=direct

numproc visualization npv1 -scalarfunction=u -subdivision=2 -nolineartexture
