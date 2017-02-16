from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.utils import NgsPickler, NgsUnpickler
import io


def test_pickle_gridfunction_real():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh,order=3,dirichlet=[1,2,3,4])
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v))

    f = LinearForm(fes)
    f += SymbolicLFI(1*v)

    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    with io.BytesIO() as f:
        pickler = NgsPickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = NgsUnpickler(f)
        u2 = unpickler.load()

    assert sqrt(Integrate((u-u2)*(u-u2),mesh)) < 1e-14


def test_pickle_gridfunction_complex():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh,order=3,complex=True,dirichlet=[1,2,3,4])
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v))

    f = LinearForm(fes)
    f += SymbolicLFI(1j*v)

    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    with io.BytesIO() as f:
        pickler = NgsPickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = NgsUnpickler(f)
        u2 = unpickler.load()
    error = sqrt(Integrate(Conj(u-u2)*(u-u2),mesh))
    assert error.real < 1e-14 and error.imag < 1e-14
