from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *
from math import pi

from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

ngsglobals.msg_level = 1

mesh = Mesh (unit_square.GenerateMesh(maxh=0.3))

v = FESpace(type="hdivho",mesh=mesh, order=1, dirichlet=[1,2,3])

u = GridFunction (space=v)

u.vec[37]=1.0

vtk = VTKOutput(ma=mesh,coefs=[u],names=["sol"],filename="vtkout_",subdivision=3)
vtk.Do()

Draw (u)
