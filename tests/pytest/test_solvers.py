from netgen.geom2d import unit_square
from ngsolve import *
import pytest
from ngsolve.krylovspace import *

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


def test_krylovspace_solvers():
    solvers = [CGSolver, GMResSolver, MinResSolver, QMRSolver]
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes = H1(mesh, order=4, dirichlet=".*")
    u,v = fes.TnT()
    f = LinearForm(32 * (y*(1-y)+x*(1-x)) * v * dx).Assemble()
    a = BilinearForm(grad(u)*grad(v)*dx)
    c = Preconditioner(a, type="bddc")
    a.Assemble()
    u = GridFunction(fes)
    exact = 16*x*(1-x)*y*(1-y)
    for solver in solvers:
        inv = solver(mat=a.mat, pre=c)
        u.vec.data = inv * f.vec
        error = sqrt(Integrate((u-exact)*(u-exact), mesh))
        print(solver.name, ": iterations = ", inv.iterations)
        print("Error = ", error)
        assert inv.iterations < 40
        # p4 should be exact
        assert error < 1e-12



if __name__ == "__main__":
    # test_arnoldi()
    test_krylovspace_solvers()
