
from netgen.geom2d import unit_square
from ngsolve import *


def Spaces():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    spaces = { \
        "h1ho" : H1(mesh,order=3), \
        "hcurlho" : HCurl(mesh,order=3) \
    }
    return spaces

timed_objects = Spaces()
