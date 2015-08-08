from ngsolve.fem import *
from ngsolve.utils import *
mesh = Mesh("square.vol")
# vtk output for all coefs in coefs=[...] which are labeled as names=[...].
# result is written to vtkout_0, vtkout_1, ... (each .Do()-call increases the number)
vtk = VTKOutput(ma=mesh,coefs=[x*x,y+4],names=["x_square","y_plus_4"],filename="vtkout_",subdivision=2)
vtk.Do()
