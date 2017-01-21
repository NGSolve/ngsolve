
from ngsolve import *
from netgen.geom2d import unit_square

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

h1 = H1(mesh=mesh,order=3)
hcurl = HCurl(mesh=mesh,order=3)
with TaskManager():
    h1timing = h1.__timing__()
    hcurltiming = hcurl.__timing__()
print("h1: ", h1timing)
print("hcurl: ", hcurltiming)

