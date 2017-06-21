from netgen.geom2d import unit_square
from ngsolve import *
import pytest

def test_arnoldi():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes1 = H1(mesh,order=5,complex=True,dirichlet="top|bottom|left|right")
    fes2 = H1(mesh,order=2,complex=True,dirichlet="top|bottom|left|right")
    fes = FESpace([fes1,fes2])

    u1,u2 = fes.TrialFunction()
    v1,v2 = fes.TestFunction()

    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u2)*grad(v1)+u1*v2)

    m = BilinearForm(fes)
    m += SymbolicBFI(u1*v1+u2*v2)

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
        assert error < 1e-10

    Draw(laplace(u.components[0]),mesh,"laplace")



if __name__ == "__main__":
    test_arnoldi()
