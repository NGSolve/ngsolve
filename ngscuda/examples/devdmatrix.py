# file: devdmatrix.py
# date: 28.09.2022
#
# testing basic functions for dense matrix on device

from ngsolve.bla import *
from ngsolve.la import *
from ngsolve.ngscuda import *

m = 3
n = 4
k = 5

A = MatrixD(m, k)
A[:] = 1

B = MatrixD(k, n)
B[:] = 2

C = MatrixD(m, n)
C[:] = 3

Adev = CreateDevMatrix(A)
Bdev = CreateDevMatrix(B)
Cdev = CreateDevMatrix(C)


# scale matrix
print("Scale Matrix")
Cdev.Scale(3)
# print(Cdev)

# add two matrices
print("Add Matrix")
Bdev.Add(Bdev)
# print(Bdev)

# set the matrix zero
print("Set Matrix to zero")
Cdev.SetZero()
# print(Cdev)

# multiply matrix-vector
print("Matrix-Vector Multiplication")
x = Adev.CreateRowVector()
x.FV()[:] = 1

y = Adev.CreateColVector()

x.UpdateDevice()
y.UpdateDevice()
Adev.Mult(x, y)
# print(y)
