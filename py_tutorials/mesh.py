import sys
sys.path.append("/opt/netgen/lib")

from libngspy import *


# pde = ngsolve.PDE("../pde_tutorial/d1_square.pde")
# mesh = pde.Mesh()

# mesh = ngcomp.Mesh("square.vol")
mesh = ngcomp.Mesh("cube.vol.gz")

# from tkinter import filedialog
# filename = filedialog.askopenfilename(filetypes=[("vol-files","*.vol *.vol.gz")])
# mesh = ngcomp.Mesh(filename)


for i in mesh.Elements(ngcomp.VOL):
    print (i)
    print (mesh[i])
    print (mesh[i].vertices)


