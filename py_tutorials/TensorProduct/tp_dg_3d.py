from TensorProductTools import *
from make_seg_mesh import *
from ngsolve.comp import *
from ngsolve import *
from ngsolve.comp import TensorProductFESpace, Transfer2StdMesh

mesh1 = Mesh(SegMesh(32,0,1))
mesh2 = Mesh(SegMesh(40,0,1))

tpmesh = Mesh(MakeTensorProductMesh(mesh1,mesh2))
Draw(tpmesh)

n=4
m=4
ngsglobals.numthreads = 24

fesx = L2(mesh1,order=n,flags={'dgjumps':True})
fesy = L2(mesh2,order=m,flags={'dgjumps':True})
tpfes = TensorProductFESpace([fesx,fesy],{'dgjumps':True})

fes = L2(tpmesh,order=n)

u = tpfes.TrialFunction()
v = tpfes.TestFunction()

vx = v.Operator("gradx")
vy = v.Operator("grady")

b = (-1)*CoefficientFunction( (x1-0.5,-x+0.5) )


uin = IfPos((1.0-x1)*(x1-0.5),sin((x1-0.5)*2*3.14159),0) + IfPos((x1-0.0)*(0.5-x1),sin((x1)*2*3.14159),0)


gradv = CoefficientFunction((vx,vy))

a = BilinearForm(tpfes)

n = specialcf.normal(tpmesh.dim)
bn = b*n

a += SymbolicBFI ( -u * b*gradv )
a += SymbolicBFI ( (bn) *IfPos(bn, u, u.Other( )) * (v-v.Other()), VOL, skeleton=True)
a += SymbolicBFI ( (bn) *IfPos(bn, u, u.Other(bnd = CoefficientFunction(uin) )) * (v), BND, skeleton=True)


u = GridFunction(tpfes)
v = GridFunction(tpfes)

uu = GridFunction(fes)

u.Set(exp(-200*(x-0.3)*(x-0.3)-200*(x1-0.3)*(x1-0.3)))

Transfer2StdMesh(u,uu)
Draw(uu,sd=4,autoscale=False)

h = u.vec.CreateVector()
print('To start the simulation type Run(n_steps)!')

def Step():
    a.Apply(u.vec,v.vec)
    h.data = 0.000625*v.vec.data
    tpfes.SolveM(rho=CoefficientFunction(1), vec = h)
    u.vec.data-=h.data
    Transfer2StdMesh(u,uu)
    Redraw()

def Run(nsteps):
    with TaskManager():
        for i in range(nsteps):
            print("Step ",i+1, "/",nsteps)
            Step()
Run(100)
for t in Timers():
    print(t["counts"], t["time"], t["name"])