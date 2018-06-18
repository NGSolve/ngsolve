from netgen import meshing
from netgen.meshing import Element0D, Element1D, MeshPoint, Pnt
from ngsolve import *

def test_hdg1d():
    m = meshing.Mesh(dim=1)
    nel = 10
    pnums = []
    for i in range(0, nel+1):
        pnums.append (m.Add (MeshPoint (Pnt(i/nel, 0, 0))))

    for i in range(0,nel):
        m.Add (Element1D ([pnums[i],pnums[i+1]], index=1))

    m.Add (Element0D (pnums[0], index=1))
    m.Add (Element0D (pnums[nel], index=1))

    mesh = Mesh(m)
    order = 2
    V1 = L2(mesh, order=order)
    V2 = FacetFESpace(mesh, order=0, dirichlet=[1])
    V = FESpace([V1,V2])

    u, uhat = V.TrialFunction()
    v, vhat = V.TestFunction()

    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size

    a = BilinearForm(V)
    a += SymbolicBFI(grad(u)*grad(v))
    a += SymbolicBFI(-grad(u)*n*(v-vhat),element_boundary=True)
    a += SymbolicBFI(-grad(v)*n*(u-uhat),element_boundary=True)
    alpha = 4
    a += SymbolicBFI(alpha*order*order/h * (u-uhat)*(v-vhat),element_boundary=True)

    f = LinearForm(V)
    f += SymbolicLFI(1*v)

    f.Assemble()
    a.Assemble()
    u = GridFunction(V)
    u.vec.data = a.mat.Inverse(V.FreeDofs(),inverse="sparsecholesky")*f.vec
    
    exsol = 0.5*x*(1-x)
    l2error = Integrate(InnerProduct(u.components[0]-exsol,u.components[0]-exsol),mesh)
    assert l2error < 1e-15

if __name__ == "__main__":
    test_hdg1d()
