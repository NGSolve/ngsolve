from ngsolve import *
import netgen.geom2d

mesh = Mesh (netgen.geom2d.unit_square.GenerateMesh(maxh=0.1))

v = FESpace ("h1ho", mesh, order=4, dirichlet=[1])
u = GridFunction (v)

data = [u.vec, u.vec, u.vec]

import pickle
# store same vector three times
pickler = pickle.Pickler(open ("1.dat", "wb"))
pickler.dump ([data])
del pickler

# NgsPickler recognizes redundancies -> smaller files
pickler = NgsPickler(open ("2.dat", "wb"))
pickler.dump ([data])
del pickler
