# file: devdmatrix.py
# date: 20.09.2022
#
# testing basic functions for dense matrix on device

from ngsolve.bla import *
from ngsolve.la import *
from ngsolve.ngscuda import *

m = 4 
n = 5
k = 7

A = MatrixD(m, n)
B = MatrixD(n, k)
C = MatrixD(m, k)
for i in range(m):
    for j in range(n):
        A[i,j] = i+2*j

for i in range(n):
    for j in range(k):
        B[i,j] = 2+2*i+j

for i in range(m):
    for j in range(k):
        C[i,j] = 1

Adev = CreateDevMatrix(A)
Bdev = CreateDevMatrix(B)
Cdev = CreateDevMatrix(C)

# set the matrix zero
Cdev.SetZero()

# scale matrix
Adev.Scale(3)

# add two matrices
Adev.Add(Adev)

# multiply two matrices
# WARNING: results of mult still wrong. will be fixed soon
Adev.Mult(Bdev, Cdev)
