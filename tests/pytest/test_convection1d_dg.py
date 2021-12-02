# generate a 1D mesh
from netgen import meshing
from netgen.csg import Pnt
from ngsolve import *
    
def test_convection1d_dg():
    m = meshing.Mesh()
    m.dim = 1
    nel = 20
    pnums = []
    for i in range(0, nel+1):
        pnums.append (m.Add (meshing.MeshPoint (Pnt(i/nel, 0, 0))))

    for i in range(0,nel):
        m.Add (meshing.Element1D ([pnums[i],pnums[i+1]], index=1))

    m.Add (meshing.Element0D (pnums[0], index=1))
    m.Add (meshing.Element0D (pnums[nel], index=2))
    m.AddPointIdentification(pnums[0],pnums[nel],identnr=1,type=meshing.IdentificationType.PERIODIC)

    mesh = Mesh (m)

    fes = L2(mesh, order=4)

    u = fes.TrialFunction()
    v = fes.TestFunction()

    b = CoefficientFunction(1)
    bn = b*specialcf.normal(1)

    a = BilinearForm(fes)
    a += SymbolicBFI (-u * b*grad(v))
    a += SymbolicBFI (bn*IfPos(bn, u, u.Other()) * v, element_boundary=True)

    u = GridFunction(fes)
    pos = 0.5
    u0 = exp (-100 * (x-pos)*(x-pos) )
    u.Set(u0)

    w = u.vec.CreateVector()

    t = 0
    tau = 1e-4
    tend = 1

    with TaskManager():
        while t < tend-tau/2:
            a.Apply (u.vec, w)
            fes.SolveM (rho=CoefficientFunction(1), vec=w)
            u.vec.data -= tau * w
            t += tau
            
    l2error = sqrt(Integrate((u-u0)*(u-u0),mesh))
    print(l2error)
    assert l2error < 1e-2
