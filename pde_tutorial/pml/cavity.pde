geometry = cavity.in2d
mesh = cavity.vol.gz

define constant one = 1
define constant pml_r = 0.8
# define constant pml_xmin = -0.8
# define constant pml_xmax = 0.8
# define constant pml_ymin = -1
# define constant pml_ymax =  0.8
define constant pml_alpha = 2

define constant geometryorder = 5
define constant hpref = 4

define fespace v -type=h1ho -order=5 -complex  -dirichlet=[3]

define gridfunction u -fespace=v  -multidim=100

define bilinearform a -fespace=v -symmetric
laplace one  -definedon=[1]
PML_laplace one -definedon=[2]

define bilinearform m -fespace=v -symmetric
mass one  -definedon=[1]
PML_mass one -definedon=[2]


numproc evp np2 -bilinearforma=a -bilinearformm=m -gridfunction=u -num=200  -shift=400 -filename=eigen.out
