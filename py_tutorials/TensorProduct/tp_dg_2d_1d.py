from TensorProductTools import *
from make_seg_mesh import *
from ngsolve.comp import *
from ngsolve import *
from ngsolve.comp import TensorProductFESpace, Transfer2StdMesh, SymbolicTPBFI
from netgen.geom2d import unit_square


mesh1 = Mesh(unit_square.GenerateMesh(maxh=0.15))
mesh2 = Mesh(SegMesh(20,0,1,1))

tpmesh = Mesh(MakeTensorProductMesh(mesh1,mesh2))
Draw(tpmesh)

n=3
m=3
ngsglobals.numthreads = 24

fesx = L2(mesh1,order=n,flags={'dgjumps':True})
fesy = L2(mesh2,order=m,flags={'dgjumps':True})
tpfes = TensorProductFESpace([fesx,fesy],{'dgjumps':True})

fes = L2(tpmesh,order=n)

u = tpfes.TrialFunction()
v = tpfes.TestFunction()

vx = v.Operator("gradx")
vy = v.Operator("grady")

b = CoefficientFunction( (0,0,1) )

uin = CoefficientFunction(0.0)

gradv = CoefficientFunction((vx,vy))

a = BilinearForm(tpfes)

n = specialcf.normal(tpmesh.dim)
bn = b*n

a += SymbolicTPBFI ( -u * b*gradv )
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other(bnd = uin )) * (v-v.Other(bnd = 0.0)), VOL, skeleton=True)
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other(bnd = uin )) * (v), BND, skeleton=True)


u = GridFunction(tpfes)
v = GridFunction(tpfes)


uu = GridFunction(fes)

u.Set(exp(-70*(x-0.125)*(x-0.125)-70*(y-0.125)*(y-0.125)-70*(x1-0.75)*(x1-0.75)))
Transfer2StdMesh(u,uu)
Draw(uu,sd=3,autoscale=False)

h = u.vec.CreateVector()
print('To start the simulation type Run(n_steps)!')
def Step():
    a.Apply(u.vec,v.vec)
    h.data = 0.001*v.vec.data
    tpfes.SolveM(rho=CoefficientFunction(1), vec = h)
    u.vec.data-=h.data
    Transfer2StdMesh(u,uu)
    Redraw()

def Run(nsteps):
    with TaskManager():
        for i in range(nsteps):
            print("Step ",i+1, "/",nsteps)
            Step()

Run(10000)
for t in Timers():
    print(t["counts"], t["time"], t["name"])