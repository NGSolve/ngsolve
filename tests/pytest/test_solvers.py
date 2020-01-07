from netgen.geom2d import unit_square
from ngsolve import *
import pytest

def test_arnoldi():
    SetHeapSize (10*1000*1000)
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.25))
    fes1 = L2(mesh,order=8,complex=True,dirichlet="top|bottom|left|right")
    fes2 = H1(mesh,order=10,complex=True,dirichlet="top|bottom|left|right")
    fes = FESpace([fes1,fes2])

    u1,u2 = fes.TrialFunction()
    v1,v2 = fes.TestFunction()

    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u2)*grad(v2)+u1*v1)

    m = BilinearForm(fes)
    m += SymbolicBFI(u1*v2+u2*v1)

    u = GridFunction(fes,multidim=40)

    a.Assemble()
    m.Assemble()

    lam = ArnoldiSolver(a.mat,m.mat,fes.FreeDofs(),u.vecs,1)
    print("ev = ", lam)
    Draw(u.components[0])
    Draw(u.components[1])

    evec = GridFunction(fes2)

    def laplace(gf):
        hesse = gf.Operator("hesse")
        return hesse[0]+hesse[3]

    for i in range(5):
        evec.vec.data = u.components[1].vecs[i]
        error = Integrate(Norm(laplace(evec) + lam[i].real*lam[i].real*evec),mesh)
        print("error[",i,"] = ",error)
        assert error < 1e-7

    Draw(laplace(evec),mesh,"laplace")

def test_newton_with_dirichlet():
    mesh = Mesh (unit_square.GenerateMesh(maxh=0.3))
    V = H1(mesh, order=3, dirichlet=[1,2,3,4])
    u,v = V.TnT()
    a = BilinearForm(V)
    a += (grad(u) * grad(v) + 3*u**3*v- 1 * v)*dx
    gfu = GridFunction(V)
    dirichlet = GridFunction(V)
    dirichlet.Set(0)
    newton = solvers.Newton(a, gfu, dirichletvalues=dirichlet.vec)


if __name__ == "__main__":
    test_arnoldi()
