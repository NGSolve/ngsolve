import pytest
from ngsolve import *

def getArray():
    a = ArrayI(10)
    for i in range(10):
        a[i] = i
    return a

def test_array_iterator():
    i = 0
    for n in getArray():
        i += 1
        assert n == i
