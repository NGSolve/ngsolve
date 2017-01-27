import netgen.meshing as meshing
from netgen.csg import *
from math import pi
from ngsolve import *
# generate a 1D mesh

def test_1dlaplace():
    m = meshing.Mesh()
    m.dim = 1

    nel = 20

    pnums = []
    for i in range(0, nel+1):
        pnums.append (m.Add (meshing.MeshPoint (meshing.Pnt(i/nel, 0, 0))))

    for i in range(0,nel):
        m.Add (meshing.Element1D ([pnums[i],pnums[i+1]], index=1))

    m.Add (meshing.Element0D (pnums[0], index=1))
    m.Add (meshing.Element0D (pnums[nel], index=2))


    mesh = Mesh(m)

    order=10

    fes = H1(mesh, order=order,dirichlet=[1,2])

    u = fes.TrialFunction()
    v = fes.TestFunction()

    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    ugf = GridFunction(fes)
    uex = GridFunction(fes)
    uex.Set(sin(pi*x))

    a = BilinearForm(fes)
    a += SymbolicBFI ( grad(u) * grad(v) )

    l = LinearForm(fes)
    l += SymbolicLFI ( pi*pi*uex*v )


    l.Assemble()
    a.Assemble()
    invmat = a.mat.Inverse(fes.FreeDofs())
    ugf.vec.data = invmat*l.vec
    error = sqrt(Integrate( (ugf-uex)*(ugf-uex) ,mesh ))
    assert error <= 1e-14
