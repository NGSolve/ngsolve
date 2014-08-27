import sys
sys.path.append("/opt/netgen/lib")

from libngspy import *
from ngstd import *
from ngbla import *
from ngfem import *
from ngcomp import *
from ngsolve import *


pde = PDE()
pde.Load ("../pde_tutorial/d1_square.pde")
mesh = pde.Mesh()

# mesh = Mesh("square.vol")
# mesh = Mesh("cube.vol.gz")

# from tkinter import filedialog
# filename = filedialog.askopenfilename(filetypes=[("vol-files","*.vol *.vol.gz")])
# mesh = Mesh(filename)


for i in mesh.Elements(VOL):
    print (i)
    print (mesh[i])
    print (mesh[i].vertices)


