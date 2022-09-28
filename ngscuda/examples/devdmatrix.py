# file: devdmatrix.py
# date: 20.09.2022
#
# testing basic functions for dense matrix on device

from ngsolve.bla import *
from ngsolve.la import *
from ngsolve.ngscuda import *

m = 3 
n = 4
k = 5

x = UnifiedVector(k * [1])
y = UnifiedVector(m)

A = MatrixD(m, k)
B = MatrixD(k, n)
C = MatrixD(m, n)
for i in range(m):
    for j in range(k):
        A[i,j] = i+2*j

for i in range(k):
    for j in range(n):
        B[i,j] = 2+2*i+j

for i in range(m):
    for j in range(n):
        C[i,j] = 1

print(A)
print(B)
print(C)

Adev = CreateDevMatrix(A)
Bdev = CreateDevMatrix(B)
Cdev = CreateDevMatrix(C)

# set the matrix zero
print("set zero")
Cdev.SetZero()
print(Cdev)

# scale matrix
print("scale")
Adev.Scale(3)
print(Adev)

# add two matrices
print("add")
Adev.Add(Adev)
print(Adev)

# multiply matrix-vector
print("mult")
x.UpdateDevice()
y.UpdateDevice()
Adev.Mult(x, y)
# y = Adev * x
print(y)

print("multadd")
Adev.MultAdd(3, x, y)
print(y)
