from netgen.csg import *

geo = CSGeometry()
left = Plane(Pnt(0,0,0),Vec(-1,0,0))
right = Plane(Pnt(1,0,0),Vec(1,0,0))
brick = OrthoBrick(Pnt(-1,0,0),Pnt(2,1,1)).bc("outer")
cube = brick * left * right
geo.Add(cube)
geo.PeriodicSurfaces(left,right)

from ngsolve import *
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

fes = Periodic(H1(mesh,order=3,dirichlet="outer"))

u,v = fes.TrialFunction(), fes.TestFunction()

a = BilinearForm(fes,symmetric=True)
a += SymbolicBFI(grad(u) * grad(v))

f = LinearForm(fes)
f += SymbolicLFI(x*v)

u = GridFunction(fes,"u")

with TaskManager():
    a.Assemble()
    f.Assemble()
    u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(u)

from ngsolve.internal import *
viewoptions.clipping.enable = 1
visoptions.clipsolution = "scal"
