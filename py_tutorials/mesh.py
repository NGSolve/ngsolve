# import sys
# import os
# from os import environ
# sys.path.append(environ['NETGENDIR']+"/../lib")

# from libngspy import *
from ngsolve import *


# pde = solve.PDE("../pde_tutorial/d1_square.pde")
# mesh = pde.Mesh()

# mesh = comp.Mesh("square.vol")
# mesh = comp.Mesh("../pde_tutorial/cube.vol")

from tkinter import filedialog
filename = filedialog.askopenfilename(filetypes=[("vol-files","*.vol *.vol.gz")])
mesh = comp.Mesh(filename)


for i in mesh.Elements(comp.VOL):
    print (i)
    print (mesh[i])
    print (mesh[i].vertices)


