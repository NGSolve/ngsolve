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


mesh = Mesh("square.vol")

v = FESpace ("h1ho", mesh, { "order" : 3, "dirichlet" : [1] })
v.Update()

u = GridFunction (space=v)
u.Update()

f = LinearForm (v,"lff", { } )
f.Add (CreateLFI ("source", 2, ConstantCF(1)))
f.Assemble(1000000)

a = BilinearForm (v,"bfa", { } )
a.Add (CreateBFI ("mass", 2, ConstantCF(1)))
a.Add (CreateBFI ("laplace", 2, ConstantCF(1)))
a.Assemble(1000000)

inv = a.mat.Inverse(v.FreeDofs())
u.vec.data = inv * f.vec





