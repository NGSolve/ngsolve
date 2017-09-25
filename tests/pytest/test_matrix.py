import pytest
from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import numpy as np

def test_matrix():
    n = 10
    m = 15
    a = Matrix(n,m)
    assert a.Height() == n
    assert a.h == n
    assert a.Width() == m
    assert a.w == m

    a[:] = 99
    a[5,:] = 55
    a[2,4] = 1
    assert a[1,2] == 99
    assert a[2,4] == 1
    assert a[5,4] == 55
    a[:,3] = 33
    assert a[2,3] == 33

def test_matrix_numpy():
    try:
        import numpy
    except:
        pytest.skip("could not import numpy")
    n = 10
    m = 15
    a = Matrix(n,m)
    assert a.Height() == n
    assert a.h == n
    assert a.Width() == m
    assert a.w == m

    a[:] = 99
    a[2,4] = 1
    assert a[1,2] == 99
    assert a[2,4] == 1

    b = a.NumPy()
    assert b[1,2] == 99
    b[1,2] = 40
    assert a[1,2] == 40
    a[1,2] = 20
    assert b[1,2] == 20

    c = Matrix(n, m, True)
    c[0,1] = 1j
    d = c.NumPy()
    assert d[0,1] == c[0,1]
    d[0,1] = 1+3j
    assert d[0,1] == c[0,1]

def test_sparsematrix_access():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    fes = H1(mesh,dim=2)
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(InnerProduct(u,v))
    a.Assemble()
    assert np.linalg.norm(np.array(a.mat[1,1]) - np.array([[0.00550458,0],[0,0.00550458]])) < 1e-8

    fes = HCurl(mesh,dim=2,complex=True)
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(InnerProduct(u,v))
    a.Assemble()
    assert abs(a.mat[1,1][0,0] - (0.33333333333325+0j)) < 1e-8

    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    fes = H1(mesh,dim=3)
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(InnerProduct(u,v))
    a.Assemble()
    assert np.linalg.norm(np.array(a.mat[1,1]) - np.array([[0.0002156,0.,0.],[0.,0.0002156,0.],[0.,0.,0.0002156]])) < 1e-8

    fes = HCurl(mesh,dim=3,complex=True)
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(InnerProduct(u,v))
    a.Assemble()
    assert abs(a.mat[1,1][0,0] - (0.016666666666666604+0j)) < 1e-8

