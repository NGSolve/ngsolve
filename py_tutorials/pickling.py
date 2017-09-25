from ngsolve import *
import netgen.geom2d

mesh = Mesh (netgen.geom2d.unit_square.GenerateMesh(maxh=0.1))

v = FESpace ("h1ho", mesh, order=4, dirichlet=[1])
v2 = L2(mesh,order=2)
u = GridFunction (v)
u2 = GridFunction(v)
vec = u.vec
data = [v,v2,u,u2,u.vec]

import pickle
pickler = pickle.Pickler(open ("1.dat", "wb"))
pickler.dump (data)
del pickler

unpickler = pickle.Unpickler(open("1.dat","rb"))
fes,fes2,w,w2,z = unpickler.load()

assert fes.mesh is fes2.mesh
assert w.space is w2.space

assert len(z) == len(u.vec)
for i in range(len(u.vec)):
    assert u.vec[i] == z[i]
