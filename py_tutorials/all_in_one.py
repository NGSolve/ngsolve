import sys
import os
from os import environ
sys.path.append(environ['NETGENDIR']+"/../lib")

from libngspy import *

# being sloppy ....
from libngspy.ngstd import *
from libngspy.ngbla import *
from libngspy.ngfem import *
from libngspy.ngcomp import *
from libngspy.ngsolve import *


import numpy as np

mesh = Mesh("square.vol")

v = FESpace ("h1ho", mesh, { "order" : 6, "dirichlet" : [1] })
v.Update()

u = GridFunction (v)
u.Update()

f = LinearForm (v)
f.Add (LFI ("source", 2, ConstantCF(1)))
f.Assemble()

a = BilinearForm (v, flags = { "symmetric" : True })
a.Add (BFI ("mass", 2, ConstantCF(1)))
a.Add (BFI ("laplace", 2, ConstantCF(1)))
a.Assemble()

inv = a.mat.Inverse(v.FreeDofs())
u.vec.data = inv * f.vec

sampling = [ (x,y,u(x,y)) for x in np.linspace(0,1,6) for y in np.linspace(0,1,6)]



