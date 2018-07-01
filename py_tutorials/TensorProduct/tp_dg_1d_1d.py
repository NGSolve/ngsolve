from ngsolve.TensorProductTools import *
from ngsolve.comp import *
from ngsolve import *
import netgen.gui

mesh1 = Mesh(SegMesh(20,0,1,periodic=True) )
mesh2 = Mesh(SegMesh(20,0,1,periodic=True) )

tpmesh = Mesh(MakeTensorProductMesh(mesh1,mesh2))
Draw(tpmesh)

n=5
m=5

fesx = L2(mesh1,order=n)
fesy = L2(mesh2,order=m)

tpfes = TensorProductFESpace([fesx,fesy])

fes = L2(tpmesh,order=n)
u = tpfes.TrialFunction()
v = tpfes.TestFunction()
vx = v.Operator("gradx")
vy = v.Operator("grady")
b = (-1)*CoefficientFunction( (0.25,0.5) )


uin = ProlongateCoefficientFunction(IfPos((1.0-x)*(x-0.5),sin((x-0.5)*2*3.14159),0) + IfPos((x-0.0)*(0.5-x),sin((x)*2*3.14159),0) , 0 , tpfes)

gradv = CoefficientFunction((vx,vy))

a = BilinearForm(tpfes)

n = CoefficientFunction((ProlongateCoefficientFunction(specialcf.normal(1),1,tpfes),ProlongateCoefficientFunction(specialcf.normal(1),0,tpfes)))
bn = b*n

a += SymbolicTPBFI ( (-u * b*gradv).Compile() )
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other( )) * (v-v.Other()).Compile(), VOL, skeleton=True)
a += SymbolicTPBFI ( (bn) *IfPos(bn, u, u.Other(bnd = CoefficientFunction(uin) )) * (v).Compile(), BND, skeleton=True)



u = GridFunction(tpfes)
v = GridFunction(tpfes)

uu = GridFunction(fes)
print("Setting")
u.Set(exp( ProlongateCoefficientFunction(-200*(x-0.4)*(x-0.4), 1, tpfes)+ ProlongateCoefficientFunction(-200*(x-0.4)*(x-0.4),0,tpfes) ) )
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
Run(100)
for t in Timers():
    print(t["counts"], t["time"], t["name"])