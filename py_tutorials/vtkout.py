from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *
from math import pi

from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

ngsglobals.msg_level = 4

mesh = Mesh (unit_square.GenerateMesh(maxh=0.3))

v = FESpace(type="hdivho",mesh=mesh, order=1, dirichlet=[1,2,3])

u = GridFunction (space=v)

vtk = VTKOutput(mesh=mesh,coefs=[u],names=["sol"],filename="vtkout",subdivision=3,floatsize="single",legacy=False)
for i in range(13,26):
    u.vec[:]=0.0
    u.vec[i]=1.0
    outputfilename = vtk.Do(time=i/100)
    print("VTKOutput has been generated in", outputfilename)

Draw (u)
