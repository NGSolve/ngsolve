from TensorProductTools import *
from make_seg_mesh import *
from ngsolve.comp import *
from ngsolve import *
from ngsolve.comp import TensorProductFESpace, Transfer2StdMesh, SymbolicTPBFI

mesh1 = Mesh(SegMesh(20,0,1,1))
mesh2 = Mesh(SegMesh(20,0,1,1))

tpmesh = Mesh(MakeTensorProductMesh(mesh1,mesh2))
Draw(tpmesh)

n=5
m=5

fesx = L2(mesh1,order=n,flags={'dgjumps':True})
fesy = L2(mesh2,order=m,flags={'dgjumps':True})

tpfes = TensorProductFESpace([fesx,fesy],{'dgjumps':True})

fes = L2(tpmesh,order=n)
u = tpfes.TrialFunction()
v = tpfes.TestFunction()
vx = v.Operator("gradx")
vy = v.Operator("grady")
#b = (-1)*CoefficientFunction( (x1-0.5,-x+0.5) )
b = (-1)*CoefficientFunction( (0.25,0.5) )
#b = CoefficientFunction( (0.3,1) )


uin = IfPos((1.0-x1)*(x1-0.5),sin((x1-0.5)*2*3.14159),0) + IfPos((x1-0.0)*(0.5-x1),sin((x1)*2*3.14159),0)
#uin = CoefficientFunction(0.4*sin(x1*3.14159) + 0.4*sin(x*3.14159) )

gradv = CoefficientFunction((vx,vy))

a = BilinearForm(tpfes)

n = specialcf.normal(tpmesh.dim)
bn = b*n

a += SymbolicTPBFI ( (-u * b*gradv).Compile() )
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other( )) * (v-v.Other()).Compile(), VOL, skeleton=True)
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other(bnd = CoefficientFunction(uin) )) * (v).Compile(), BND, skeleton=True)



u = GridFunction(tpfes)
v = GridFunction(tpfes)

uu = GridFunction(fes)
print("Setting")
u.Set(exp(-200*(x-0.4)*(x-0.4)-200*(x1-0.4)*(x1-0.4)))
print("Done")
Transfer2StdMesh(u,uu)

Draw(uu,sd=2,autoscale=False)

h = u.vec.CreateVector()
#u.vec[:] = 1.0
print('To start the simulation type Run(n_steps)!')

def Step():
    a.Apply(u.vec,v.vec)
    h.data = 0.00125*v.vec.data
    tpfes.SolveM(rho=CoefficientFunction(1), vec = h)
    u.vec.data-=h.data
    #Redraw()

def Run(nsteps):
    count = 0
    with TaskManager():
        for i in range(nsteps):
            print("Step ",i+1, "/",nsteps)
            Step()
            count += 1
            if count % 1 == 0:
                Transfer2StdMesh(u,uu)
                count = 0
                Redraw()

Transfer2StdMesh(u,uu)
Run(1500)
for t in Timers():
    print(t["counts"], t["time"], t["name"])