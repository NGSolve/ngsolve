import pytest
from ngsolve import *
from netgen.geom2d import *
from netgen.csg import unit_cube

ngsglobals.msg_level=0

def test_mass_l2():
    mesh2 = Mesh(unit_square.GenerateMesh(maxh=2))
    mesh3 = Mesh(unit_cube.GenerateMesh(maxh=2))
    meshes = [mesh2, mesh3]

    for mesh in meshes:
        l2 = L2(mesh, order=3)
        c1 = FESpace([l2,l2])
        fes = FESpace([l2,FESpace([l2,l2])])

        mass = BilinearForm(fes)
        u = fes.TrialFunction()
        v = fes.TestFunction()
        mass += SymbolicBFI(u[0]*v[0] + u[1][0]*v[1][0] + u[1][1]*v[1][1])

        mass.Assemble()
        m = mass.mat

        # check if m is diagonal
        for d in range(fes.ndof):
            m[d,d] = 0.0

        assert Norm(m.AsVector())/fes.ndof<1e-15

# check if component gf keeps gf alive
def test_component_keeps_alive():
    mesh = Mesh(unit_square.GenerateMesh(maxh=2))
    assert mesh.ne == 2
    fes1 = H1(mesh)
    fes = FESpace([fes1,fes1])
    gf = GridFunction(fes)
    gf.vec[:] = 1
    gf1, gf2 = gf.components
    gf11 = gf.components[0]
    assert gf1 == gf11 # should return the same object
    assert len(gf.vec) == 8
    assert len(gf1.vec) == 4
    mesh.Refine()
    gf.Update()
    # check if gf update calls gf1 update as well
    assert len(gf.vec) == 18
    assert len(gf1.vec) == 9
    del gf
    # check if vec is still alive
    assert len(gf1.vec) == 9
    gf1.Set(0)
    for val in gf1.vec:
        assert val == 0.0

# test partially definedon number spaces
def test_number_definedon():
    geo = CSG2d()

    inner = Rectangle((0.3, 0.3), (0.4, 0.4)).Mat("inner")
    outer = Rectangle((0, 0), (1, 1))
    geo.Add((outer-inner).Mat("outer"))
    geo.Add(inner)

    mesh = Mesh(geo.GenerateMesh())

    Draw(mesh)

    h1 = H1(mesh)
    number = NumberSpace(mesh) * NumberSpace(mesh)
    number2 = NumberSpace(mesh, definedon="inner") * NumberSpace(mesh, definedon="outer")
    u = h1.TestFunction()
    t11, t12 = number.TrialFunction()
    t21, t22 = number2.TrialFunction()

    a1 = BilinearForm(testspace=h1, trialspace=number)
    a1 += u * t11 * dx("inner")
    a1 += u * t12 * dx("outer")

    a2 = BilinearForm(testspace=h1, trialspace=number2)
    a2 += u * t21 * dx("inner")
    a2 += u * t22 * dx("outer")

    a1.Assemble()
    a2.Assemble()

    for i in range(a1.mat.height):
        for j in range(a1.mat.width):
            assert a1.mat[i,j] == pytest.approx(a2.mat[i,j])

def test_vector_interleaved():
    for mesh in [Mesh(unit_square.GenerateMesh(maxh=0.3)), Mesh(unit_cube.GenerateMesh(maxh=0.5))]:
        D = mesh.dim
        cf = CF((x*y, sin(x), z)[:D])
        makers = [lambda il: VectorH1(mesh, order=3, dirichlet="left", dirichlety="bottom", interleaved=il),
                  lambda il: VectorValued(H1(mesh, order=3, dirichlet="left"), interleaved=il)]
        for make in makers:
            fb, fi = make(False), make(True)
            nd = fi.ndof//D
            perm = [j*nd+k for k in range(nd) for j in range(D)]
            for i in range(fi.ndof):
                assert fi.CouplingType(i) == fb.CouplingType(perm[i])
                assert fi.FreeDofs()[i] == fb.FreeDofs()[perm[i]]

            sols = []
            for fes in (fb, fi):
                u,v = fes.TnT()
                a = BilinearForm(InnerProduct(Grad(u),Grad(v))*dx+u*v*dx, condense=True)
                pre = Preconditioner(a, "bddc")
                a.Assemble()
                f = LinearForm(cf*v*dx).Assemble()
                gf = GridFunction(fes)
                gf.Set(cf)
                assert Integrate((gf-cf)**2, mesh) < 1e-8
                gf.vec[:] = 0
                solvers.BVP(bf=a, lf=f, gf=gf, pre=pre, tol=1e-12, print=False)
                sols.append(gf.vec)
            for i in range(fi.ndof):
                assert sols[1][i] == pytest.approx(sols[0][perm[i]], abs=1e-10)

            with pytest.raises(Exception):
                GridFunction(fi).components[0]
            with pytest.raises(Exception):
                fi.Embedding(0)

        # dof transformations act on the blocked element-local order
        fb = VectorValued(HCurl(mesh, order=2), 2)
        fi = VectorValued(HCurl(mesh, order=2), 2, interleaved=True)
        nd = fi.ndof//2
        perm = [j*nd+k for k in range(nd) for j in range(2)]
        mcf = CF((y, x*x) + ((z,) if D==3 else ()) + (x*y, 1) + ((x,) if D==3 else ()), dims=(2,D))
        gb, gi = GridFunction(fb), GridFunction(fi)
        gb.Set(mcf); gi.Set(mcf)
        for i in range(fi.ndof):
            assert gi.vec[i] == pytest.approx(gb.vec[perm[i]], abs=1e-12)

if __name__ == "__main__":
    test_component_keeps_alive()
