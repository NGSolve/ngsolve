import pytest
from ngsolve import *

class MyMatrix(BaseMatrix):
    def IsComplex(self):
        return False

    # y = self * x
    def Mult(self, x, y):
        for i in range(len(x)):
            y[i] = sum(x[:i])

    # y += alpha * self * x
    def MultTransAdd(self, s, x, y):
        for i in range(len(x)):
            y[i] = s*sum(x[i:])

    def CreateColVector(self):
        return CreateVVector(5)

def test_derive_basematrix():
    m = MyMatrix()

    x = CreateVVector(5)
    y = CreateVVector(5)

    for i in range(len(x)):
        x[i]=1+i

    x *= 10
    for i in range(len(x)):
        assert x[i] == 10+10*i

    y.data = 2.0*m*x
    assert list(y) == [0.0, 20.0, 60.0, 120.0, 200.0]

    y.data = 2.0*m.T*x
    assert list(y) == [300.0, 280.0, 240.0, 180.0, 100.0]



def test_derive_basematrix_lifetime():
    m = MyMatrix()@MyMatrix()
    m1 = MyMatrix()

    x = CreateVVector(5)
    y = CreateVVector(5)
    y1 = CreateVVector(5)

    for i in range(len(x)):
        x[i]=1+i

    y.data = m*x
    y1.data = (m1@m1)*x
    assert list(y) == list(y1)




