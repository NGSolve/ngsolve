from netgen.csg import *
from ngsolve import *


def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(2,1,2)).bc("outer")

    core = OrthoBrick(Pnt(0,-0.05,0),Pnt(0.8,0.05,1))- \
           OrthoBrick(Pnt(0.1,-1,0.1),Pnt(0.7,1,0.9))- \
           OrthoBrick(Pnt(0.5,-1,0.4),Pnt(1,1,0.6)).maxh(0.2).mat("core")
    
    coil = (Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.3) - \
            Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.15)) * \
            OrthoBrick (Pnt(-1,-1,0.3),Pnt(1,1,0.7)).maxh(0.2).mat("coil")
    
    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)
    return geometry



ngmesh = MakeGeometry().GenerateMesh(maxh=0.5)
ngmesh.Save("coil.vol")
mesh = Mesh(ngmesh)

# curve elements for geometry approximation
mesh.Curve(5)

ngsglobals.msg_level = 5

V = HCurl(mesh, order=4, dirichlet="outer", flags = { "nograds" : True })

# u and v refer to trial and test-functions in the definition of forms below
u = V.TrialFunction()
v = V.TestFunction()

mur = { "core" : 1000, "coil" : 1, "air" : 1 }
mu0 = 1.257e-6
nu_coef = [ 1/(mu0*mur[mat]) for mat in mesh.GetMaterials() ]
print ("nu_coef=", nu_coef)

nu = DomainConstantCF (nu_coef)
a = BilinearForm(V, symmetric=True)
a += SymbolicBFI(nu*curl(u)*curl(v) + 1e-6*nu*u*v)

c = Preconditioner(a, type="bddc")
# c = Preconditioner(a, type="multigrid", flags = { "smoother" : "block" } )
a.Assemble()

f = LinearForm(V)
f += SymbolicLFI(CoefficientFunction((y,-x,0)) * v, definedon=mesh.Materials("coil"))
f.Assemble()

u = GridFunction(V)

solver = CGSolver(mat=a.mat, pre=c.mat)
u.vec.data = solver * f.vec


Draw (u.Deriv(), mesh, "B-field", draw_surf=False)
