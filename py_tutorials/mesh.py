import sys
sys.path.append("/opt/netgen/lib")

import libngspy
import ngstd
import ngbla
import ngfem
import ngcomp
import ngsolve




pde = ngsolve.PDE("../pde_tutorial/d1_square.pde")
mesh = pde.Mesh()

# mesh = Mesh("square.vol")
# mesh = Mesh("cube.vol.gz")

# from tkinter import filedialog
# filename = filedialog.askopenfilename(filetypes=[("vol-files","*.vol *.vol.gz")])
# mesh = Mesh(filename)


for i in mesh.Elements(ngcomp.VOL):
    print (i)
    print (mesh[i])
    print (mesh[i].vertices)


