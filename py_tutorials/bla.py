from ngsolve.bla import *

# number of rows
n = 10

# number of columns
m = 7

x = Vector(m)
A = Matrix(n,m)

# initialize A and x
for i in range(m):
    x[i] = i+1
    for j in range(n):
        A[j,i] = j*j+i+1
        A[j,i] = i*i+j+1


# Arithmetic
y = A*x
z = 5*A*x + 7*y - A*(x-y[0:m])


# Get third row as Vector
A[2]

# Get third row as Matrix
A[2:3]

# Get rows 4 to 7
A[4:8]

# Get columns 1,3,5 (0-based)
A[:, 1:6:2]

# Set every second element in every third column to ten
A[::2,::3] = 10

# Multiply this elements by 7
A[::2,::3] *= 7

# Reset elements in last two rows and last three columns
A[-2:,-3:] = 0

# other stuff
A[3] = x
A[:, 2] = A*x

# Complex classes
B = Matrix(10,10, complex=True)
C = Matrix(10,10, True)
w = Vector(10,True) 

for i in range(10):
    w[i] = complex(i+1, i*i)
    for k in range(10):
        B[k,i] = complex(k*k+i+1, i+1+k)

v = B*w


# Create numpy matrix from Matrix
from numpy import *
AA = asmatrix(A)
xx = asmatrix(x).T
print( AA*xx - asmatrix(A*x).T )
