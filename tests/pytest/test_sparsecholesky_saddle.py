import pytest
from ngsolve import *
from netgen.occ import *

ngsglobals.msg_level = 0

def sq(maxh): return Mesh(unit_square.GenerateMesh(maxh=maxh))
def cube(maxh): return Mesh(unit_cube.GenerateMesh(maxh=maxh))

def neumann_mean(mesh, complex=False):
    X = H1(mesh, order=2, complex=complex)*NumberSpace(mesh, complex=complex)
    (u, l), (v, m) = X.TnT()
    a = BilinearForm(X, symmetric=True)
    a += (grad(u)*grad(v) + l*v + m*u)*dx
    return a, X.FreeDofs()

def two_components():
    # two floating domains, each with its own mean value
    shape = Glue([Rectangle(1, 1).Face().Move((2*i, 0, 0)) for i in range(2)])
    shape.faces[0].name = "a"
    shape.faces[1].name = "b"
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.1))
    X = H1(mesh, order=2)*NumberSpace(mesh, definedon="a")*NumberSpace(mesh, definedon="b")
    (u, la, lb), (v, ma, mb) = X.TnT()
    a = BilinearForm(X, symmetric=True)
    a += grad(u)*grad(v)*dx
    a += (la*v + ma*u)*dx("a") + (lb*v + mb*u)*dx("b")
    return a, X.FreeDofs()

def rigid_body(mesh):
    # pure traction elasticity, rigid body motions fixed by multipliers
    dim = mesh.dim
    if dim == 2:
        rbm = [CF((1, 0)), CF((0, 1)), CF((-y, x))]
    else:
        rbm = [CF((1, 0, 0)), CF((0, 1, 0)), CF((0, 0, 1)),
               CF((-y, x, 0)), CF((0, -z, y)), CF((z, 0, -x))]
    X = VectorH1(mesh, order=2)
    for r in rbm: X = X*NumberSpace(mesh)
    u, v = X.TrialFunction(), X.TestFunction()
    a = BilinearForm(X, symmetric=True)
    a += InnerProduct(Sym(Grad(u[0])), Sym(Grad(v[0])))*dx
    for i, r in enumerate(rbm):
        a += (u[i+1]*(r*v[0]) + v[i+1]*(r*u[0]))*dx
    return a, X.FreeDofs()

def stokes(mesh, order=2, mean=False):
    V = VectorH1(mesh, order=order, dirichlet=".*" if mean else "bottom")
    Q = H1(mesh, order=order-1)
    X = V*Q*NumberSpace(mesh) if mean else V*Q
    if mean: (u, p, l), (v, q, m) = X.TnT()
    else: (u, p), (v, q) = X.TnT()
    a = BilinearForm(X, symmetric=True)
    a += (InnerProduct(grad(u), grad(v)) + div(u)*q + div(v)*p)*dx
    if mean: a += (l*q + m*p)*dx
    return a, X.FreeDofs()

cases = {
    "neumann2d": lambda: neumann_mean(sq(0.05)),
    "neumann3d": lambda: neumann_mean(cube(0.2)),
    "neumann2d_complex": lambda: neumann_mean(sq(0.05), complex=True),
    "two_components": two_components,
    "rigid_body2d": lambda: rigid_body(sq(0.05)),
    "rigid_body3d": lambda: rigid_body(cube(0.2)),
    "stokes2d": lambda: stokes(sq(0.05)),
    "stokes2d_mean": lambda: stokes(sq(0.05), mean=True),
    "stokes2d_p4": lambda: stokes(sq(0.1), order=4),
    "stokes3d": lambda: stokes(cube(0.25)),
}

def residual(a, fd, **kwargs):
    b = a.mat.CreateVector()
    b.SetRandom(1)
    proj = Projector(fd, True)
    b.data = proj*b
    x = b.CreateVector()
    x.data = a.mat.Inverse(fd, inverse="sparsecholesky", **kwargs)*b
    r = (proj*(b - a.mat*x)).Evaluate()
    return Norm(r)/Norm(b)

# default ordering, with and without the dense row detection
@pytest.mark.parametrize("flags", [{}, {"densefac": 0}])
@pytest.mark.parametrize("name", cases)
def test_saddle(name, flags):
    a, fd = cases[name]()
    a.Assemble()
    assert residual(a, fd, flags=flags) < 1e-10

def test_lastdofs():
    a, fd = neumann_mean(sq(0.05))
    a.Assemble()
    n = a.mat.height
    assert residual(a, fd, flags={"densefac": 0, "late": "off", "lastdofs": [n-1]}) < 1e-10

def test_wrong_flags():
    a, fd = neumann_mean(sq(0.2))
    a.Assemble()
    with pytest.raises(Exception):
        a.mat.Inverse(fd, inverse="sparsecholesky", flags={"late": "never"})
    with pytest.raises(Exception):
        a.mat.Inverse(fd, inverse="sparsecholesky", flags={"lastdofs": [a.mat.height]})

if __name__ == "__main__":
    for name in cases:
        a, fd = cases[name]()
        a.Assemble()
        print(f"{name:20s} ndof {sum(fd):6d}  default {residual(a, fd):.1e}"
              f"  densefac=0 {residual(a, fd, flags={'densefac': 0}):.1e}"
              f"  late=off {residual(a, fd, flags={'late': 'off'}):.1e}")
