from netgen.occ import *
from ngsolve import *


def MakeGeometry():
    box = Box((-1,-1,-1), (2,1,2))
    box.faces.name = "outer"
    core = Box((0,-0.05,0), (0.8,0.05,1)) - \
        Box((0.1,-1,0.1), (0.7,1,0.9)) - \
        Box((0.5,-1,0.4), (1,1,0.6))
    core.solids.name = "core"
    core.solids.maxh = 0.03

    coil = Cylinder((0.05,0,0.3), (0,0,1), 0.3, 0.4) - \
        Cylinder((0.05,0,0.3), (0,0,1), 0.15, 0.4)
    coil.solids.name = "coil"
    air = box - core - coil
    air.solids.name = "air"
    return Glue([core, coil, air])



mesh = MakeGeometry().GenerateMesh(maxh=0.5)

# curve elements for geometry approximation
mesh.Curve(5)

ngsglobals.msg_level = 5

fes = HCurl(mesh, order=4, nograds = True)

# u and v refer to trial and test-functions in the definition of forms below
u,v = fes.TnT()

mur = mesh.MaterialCF({ "core" : 1000 }, default=1)
mu0 = 1.257e-6
nu = 1/(mu0*mur)

# the variational equation, with the current density given in the coil
lhs = nu*curl(u)*curl(v)*dx + 1e-6*nu*u*v*dx
rhs = CF((y,0.05-x,0)) * v * dx("coil")

# Dirichlet condition on the named boundary "outer",
# solve by conjugate gradients with BDDC preconditioner
gfu = Solve(lhs==rhs, u[BND("outer")]==(0,0,0),
            preconditioners.BDDC(), solvers.CGSolver.Creator())


Draw (curl(gfu), mesh, "B-field", draw_surf=False)
