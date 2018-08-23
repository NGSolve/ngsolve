from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry
from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (1, 1), 
                 bcs = ("bottom", "right", "top", "left"))
mesh = Mesh( geo.GenerateMesh(maxh=0.2))

b = CoefficientFunction((1,2))
n = specialcf.normal(mesh.dim)
ubnd = IfPos(x,1,0)

V = L2(mesh,order=2,flags={"dgjumps": True})
a = BilinearForm(V)
u,v = V.TrialFunction(), V.TestFunction()
a += SymbolicBFI( b * grad(u) * v)
a += SymbolicBFI( IfPos( (b*n), 0, (b*n) * (u.Other()-u)) * v, element_boundary=True)
a.Assemble()

f = LinearForm(V)
f += SymbolicLFI( (b*n) * IfPos( (b*n), 0, -ubnd) * v, BND, skeleton=True)
f.Assemble()

gfu_impl = GridFunction(V)
gfu_impl.vec.data = a.mat.Inverse() * f.vec
Draw(gfu_impl,mesh,"u_implicit")
