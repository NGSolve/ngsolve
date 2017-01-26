import pytest
from ngsolve import *

def test_vector():
    n = 10
    a = Vector(n)
    assert len(a) == n
    for i in range(n):
        a[i] = i
    assert list(a) == list(range(n))
    assert list(a-a) == [0]*n

    b = Vector(n)
    b[:] = 42
    b[1:8:3] = 99

    assert list(b)== ([42,99,42]*3)+[42]

def test_vector_numpy():
    try:
        import numpy
    except:
        pytest.skip("could not import numpy")
    n = 10
    a = Vector(n)
    for i in range(len(a)):
        a[i] = i

    v = a.NumPy()
    v[0] = 234
    assert a[0] == 234
    a[1] = 100
    assert v[1] == 100

    assert list(v) == list(a)

    c = Vector(n, True)
    c[0] = 1j
    d = c.NumPy()
    assert d[0] == c[0]
    d[1] = 1+3j
    assert d[1] == c[1]
