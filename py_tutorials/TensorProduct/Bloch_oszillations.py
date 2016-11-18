from TensorProductTools import *
from make_seg_mesh import *
from ngsolve.comp import *
from ngsolve import *
from ngsolve.comp import TensorProductFESpace, Transfer2StdMesh, IntDv,IntDv2, Prolongate
from math import pi
from netgen.geom2d import unit_square
ngsglobals.numthreads = 12
ngsglobals.msg_level = 0
L = 2.0
Lk = 1
j=3
i=1.5


linear = False
if linear==False:
    epsm =  0.25*( 1-cos(x1*pi/Lk))
    depsm = 0.25*sin(x1*pi/Lk)/Lk*pi

if linear==True:
    c=1/ ( 2/(i*j) - 1/j**2 - 1/(i*j)*1/(1-i) + 2/j*1/(1-i) + i/j*1/(1-i) )
    a = 2*c*1/j
    b = -c/(j**2)
    d=i/j*1/(1-i)*c
    e=-2*i/j*1/(1-i)*c
    f=1-i/j*1/(1-i)*c

    epsmright = IfPos(1/j-x1 ,-c*x1*x1,IfPos(1/i-x1,-a*x1-b, -d*x1*x1-e*x1-f   )  )
    epsmleft =  IfPos(-1/i-x1, -d*x1*x1+e*x1-f, IfPos(-1/j-x1,  a*x1-b, -c*x1*x1)   )

    depsmright = IfPos(1/j-x1 ,-2*c*x1,IfPos(1/i-x1,-a, -2*d*x1-e   )  )
    depsmleft =  IfPos(-1/i-x1, -2*d*x1+e, IfPos(-1/j-x1,  a, -2*c*x1) )
    epsm = IfPos(x1, 0.5*epsmright, 0.5*epsmleft)
    depsm = IfPos(x1, 0.5*depsmright, 0.5*depsmleft)
def RK4(gm,am,fes,tau,depsm):
    am.Apply(gm.vec, k1m)
    h1m.data = gm.vec-tau/2.0*k1m    
    am.Apply(h1m, k2m)
    h2m.data = gm.vec-tau/2.0*k2m   
    am.Apply(h2m ,k3m)
    h3m.data = gm.vec-tau*k3m
    am.Apply(h3m,k4m)
    wm.data = k1m+2.0*k2m+2.0*k3m+k4m
    fes.SolveM(rho = CoefficientFunction(1), vec = wm)
    gm.vec.data -= tau/6.0*wm


mesh1 = Mesh(SegMesh(50,   0.0, L,0))
mesh2 = Mesh(SegMesh(50,  -Lk,  Lk,1))

tpmesh = Mesh(MakeTensorProductMesh(mesh1,mesh2))    # for visualization only, std netgen mesh!!!
Draw(tpmesh)                                         # for visualization only

orderx=3
ordery=3

# FESpaces for the problem
fesx = L2(mesh1,order=orderx,flags={'dgjumps':True})
#fesx2 = L2(mesh1,order=orderx+2,flags={'dgjumps':True})
fesy = L2(mesh2,order=ordery,flags={'dgjumps':True})
tpfes = TensorProductFESpace([fesx,fesy],{'dgjumps':True})

#tpfes2 = TensorProductFESpace([fesx2,fesy],{'dgjumps':True})
fes = L2(tpmesh,order=orderx)           # for visualization only

# Test and trial functions

u2d = tpfes.TrialFunction()
v2d = tpfes.TestFunction()

vx = v2d.Operator("gradx")
vy = v2d.Operator("grady")
gradv = CoefficientFunction((vx,vy))

gm = GridFunction(tpfes)
ggm = GridFunction(fes)                 # for visualization only

velm = GridFunction(tpfes)
velmvis = GridFunction(fes)
velm.Set(depsm)
Transfer2StdMesh(velm,velmvis)

#dV2ddx = GridFunction(tpfes)            # the gradient of the eletrco static potential
#dVdxvis = GridFunction(fes)             # and its prolongation to the to space

# Initial Distribution
x0 = 1.0
v0 = 0.0
print("Setting")
gm.Set(exp( -20*(x-x0)*(x-x0) -20*(x1-v0)*(x1-v0) ) )


# 2D transport
E = GridFunction(tpfes)
E.Set(CoefficientFunction(1.0))
bm = CoefficientFunction( (depsm,  E ) )

n = specialcf.normal(tpmesh.dim)
bnm = bm*n
am = BilinearForm(tpfes)
am += SymbolicBFI ( (-u2d * bm*gradv ))
am += SymbolicBFI ( (bnm *IfPos(bnm, u2d, u2d.Other()) * (v2d-v2d.Other())), VOL, skeleton=True)
am += SymbolicBFI ( (bnm *IfPos(bnm, u2d, u2d.Other(bnd=0.0)) * v2d), BND, skeleton=True)

Transfer2StdMesh(gm,ggm)
Draw(velmvis, tpmesh,"V_m")
Draw(ggm,tpmesh, "gm")

k1m = gm.vec.CreateVector()
k2m = gm.vec.CreateVector()
k3m = gm.vec.CreateVector()
k4m = gm.vec.CreateVector()
h1m = gm.vec.CreateVector()
h2m = gm.vec.CreateVector()
h3m = gm.vec.CreateVector()
wm  = gm.vec.CreateVector()    

#gmtotal = GridFunction(fes, multidim=50)
redrawfreq = 15
def RunRK4par(tend=2.0,tau = 1.0e-3):
    counter = 0
    mcomp = 0
    with TaskManager():
        t=0
        while t<tend:
            counter+=1
            RK4(gm,am,tpfes,tau,depsm)
            if counter % redrawfreq == 0:
                Transfer2StdMesh(gm,ggm)
                Redraw()            
            #if counter % redrawfreq == 0:
            print(t)
            t += tau

            
RunRK4par(10.1,0.00125)
for t in Timers():
    print(t["counts"], t["time"], t["name"])