from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (2, 0.41), bcs = ("wall", "outlet", "wall", "inlet") )
geo.AddCircle ( (0.2, 0.2), r=0.05, leftdomain = 0, rightdomain = 1, bc = "cyl" )
mesh = Mesh( geo.GenerateMesh(maxh = 0.08) )
order = 3
mesh.Curve(order)

V1 = HDiv ( mesh, order = order, dirichlet = "wall|cyl|inlet" )
V2 = FESpace ( "vectorfacet", mesh, order = order, dirichlet = "wall|cyl|inlet" )
Q = L2( mesh, order = order-1)
V = FESpace ([V1,V2,Q])

u, uhat, p = V.TrialFunction()
v, vhat, q = V.TestFunction()

nu = 0.001
alpha = 4

gradu = CoefficientFunction ( (grad(u),), dims=(2,2) )
gradv = CoefficientFunction ( (grad(v),), dims=(2,2) )

n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size

def tang(vec):
    return vec - (vec*n)*n

# HDG Stokes formulation:
a = BilinearForm ( V, symmetric=True)
a += SymbolicBFI ( nu*InnerProduct ( gradu, gradv) )
a += SymbolicBFI ( nu*InnerProduct ( gradu.trans * n,  tang(vhat-v) ), element_boundary=True )
a += SymbolicBFI ( nu*InnerProduct ( gradv.trans * n,  tang(uhat-u) ), element_boundary=True )
a += SymbolicBFI ( nu*alpha*order*order/h * InnerProduct ( tang(vhat-v),  tang(uhat-u) ), element_boundary=True )
a += SymbolicBFI ( -div(u)*q -div(v)*p )
a.Assemble()

m = BilinearForm(V , symmetric=True)
m += SymbolicBFI( u * v )
m.Assemble()

f = LinearForm ( V )

gfu = GridFunction(V)

# boundary values
U0 = 1.5
uin = CoefficientFunction( (U0*4*y*(0.41-y)/(0.41*0.41),0) )
gfu.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))


invstokes = a.mat.Inverse(V.FreeDofs(), inverse="umfpack")

res = f.vec.CreateVector()
res.data = f.vec - a.mat*gfu.vec
gfu.vec.data += invstokes * res

Draw( gfu.components[0], mesh, "velocity" )
Draw( Norm(gfu.components[0]), mesh, "absvel(hdiv)")

VL2 = L2(mesh, dim=mesh.dim, order=order, flags = { "dgjumps" : True } )
uL2 = VL2.TrialFunction()   # 1 x dim matrix
vL2 = VL2.TestFunction()    

gfuL2 = GridFunction(VL2)

#mixed mass matrices:
bfmixed = BilinearForm ( V, VL2, flags = { "nonassemble" : True } )
bfmixed += SymbolicBFI ( vL2*u )

bfmixedT = BilinearForm ( VL2, V, flags = { "nonassemble" : True } )
bfmixedT += SymbolicBFI ( uL2*v )

#convection operator:
vel = gfu.components[0]
convL2 = BilinearForm(VL2, flags = { "nonassemble" : True } )
convL2 += SymbolicBFI( (-InnerProduct(grad(vL2).trans * vel, uL2.trans)) )
un = InnerProduct(vel,n)
upwindL2 = IfPos(un, un*uL2, un*uL2.Other(bnd=uin))
convL2 += SymbolicBFI( InnerProduct (upwindL2, vL2-vL2.Other()), VOL, skeleton=True )
convL2 += SymbolicBFI( InnerProduct (upwindL2, vL2), BND, skeleton=True )

#convection step:
def SolveConvectionSteps(gfuvec, res, tau, steps):
    bfmixed.Apply (gfuvec, gfuL2.vec) 
    VL2.SolveM(CoefficientFunction(1), gfuL2.vec)
    conv_applied = gfuL2.vec.CreateVector()
    for i in range(steps):
        convL2.Apply(gfuL2.vec,conv_applied)
        VL2.SolveM(CoefficientFunction(1), conv_applied)
        gfuL2.vec.data -= tau/steps * conv_applied
        #Redraw()    
    bfmixedT.Apply (gfuL2.vec, res)

# initial values again:
res.data = f.vec - a.mat*gfu.vec
gfu.vec.data += invstokes * res
t = 0

tend = 1
substeps = 10
tau = 0.01

mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + tau * a.mat.AsVector()
inv = mstar.Inverse(V.FreeDofs(), inverse="umfpack")

while t < tend:
    SolveConvectionSteps(gfu.vec, res, tau, substeps)
    res.data -= mstar * gfu.vec
    gfu.vec.data += inv * res
    t += tau
    print ("t=", t)
    Redraw(blocking=True)
