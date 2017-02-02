from ngsolve import *
from netgen.geom2d import SplineGeometry
from math import sqrt

geo = SplineGeometry()

s = 0.3
Points = [(-s/2,-sqrt(3)/6*s),(s/2,-sqrt(3)/6*s),(0,sqrt(3)/3*s),
        (-s/2*4,-sqrt(3)/6*s*4),(s/2*4,-sqrt(3)/6*s*4),(0,sqrt(3)/3*s*4),
        (-s/2*6,-sqrt(3)/6*s*6),(s/2*6,-sqrt(3)/6*s*6),(0,sqrt(3)/3*s*6)]

for pnt in Points:
    geo.AddPoint(*pnt) 

for pids in [[0,1],[1,2],[2,0]]:
    geo.Append(["line"]+pids,leftdomain=0,rightdomain=1,bc="scatterer")

for pids in [[3,4],[4,5],[5,3]]:
    geo.Append(["line"]+pids,leftdomain=1,rightdomain=2)

for pids in [[6,7],[7,8],[8,6]]:
    geo.Append(["line"]+pids,leftdomain=2,rightdomain=0)


ngmesh = geo.GenerateMesh(maxh=0.04)
# ngmesh.Save("scattering.vol")
mesh = Mesh(ngmesh)
# mesh = Mesh ("scattering.vol")
p=pml.HalfSpace(point=(0,4/3*sqrt(3)*s),normal=(sqrt(3),1),alpha=0.5j)+pml.HalfSpace(point=(0,4/3*sqrt(3)*s),normal=(-sqrt(3),1),alpha=0.5j)+pml.HalfSpace(point=(0,-sqrt(3)*s*2/3),normal=(0,-1),alpha=0.5j)
print(p)
mesh.SetPML(p,definedon=2)

kx = 50
ky = 10
k = sqrt(kx*kx+ky*ky)

uin = exp (1J*kx*x+1J*ky*y)

fes = H1(mesh, complex=True, order=5, dirichlet="scatterer")
u = fes.TrialFunction()
v = fes.TestFunction()

uscat = GridFunction (fes)
uscat.Set (uin, definedon=mesh.Boundaries("scatterer"))

a = BilinearForm (fes, symmetric=True)
a += SymbolicBFI (grad(u)*grad(v) )
a += SymbolicBFI (-k*k*u*v)

f = LinearForm (fes)

a.Assemble()
f.Assemble()

res = uscat.vec.CreateVector()
res.data = f.vec - a.mat * uscat.vec
uscat.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * res


Draw (uin, mesh, "uin")
Draw (uscat, mesh, "uscat")
Draw (uin-uscat, mesh, "utot")
cf=mesh.GetPMLTrafo(2).PML_CF(2)
mesh.UnSetPML(2)
Draw (cf-Conj(cf),mesh,"scaling")




