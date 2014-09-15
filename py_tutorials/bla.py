import sys
import os
sys.path.append(os.environ["NETGENDIR"]+'/../lib')
sys.path.append(os.environ["NETGENDIR"])

from libngspy.ngbla import *

# number of rows
n = 10

# number of columns
m = 7

x = Vector(m)
A = Matrix(n,m)

#initialize A and x
for i in range(m):
    x[i] = i+1
    for j in range(n):
        A[j,i] = j*j+i+1
        A[j,i] = i*i+j+1

# Allocate result Vector
y = Vector(n)

# Write A*x to already existing Vector y,
# -> no extra memory allocation!
y.data = A*x

# Allocate new Vector with result
z = Vector(5*A*x + 7*y)

# Should be zero
print(Vector(z-12*y))

# my_expression is of type 'SumExpr'
# -> no calculation is done at this point
my_expression = x+y

# Assign x+y to z
z.data = my_expression
print(z)

