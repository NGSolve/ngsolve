import pytest
from ngsolve.ngstd import ArrayI

def getArray():
    a = ArrayI(10)
    for i in range(10):
        a[i] = i
    return a

def test_array_iterator():
    i = 0
    a = getArray()
    for n in a:
        assert n == i
        i += 1

def test_array_iterator_temporary_object():
    i = 0
    for n in getArray():
        assert n == i
        i += 1
