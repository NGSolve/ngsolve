# file: unifiedvector.py
# date: 20.09.2022
#
# testing basic functions for unifiedvectors working on device

from ngsolve.la import *
from ngsolve.ngscuda import *

n = 10

# create vectors using size
a0 = UnifiedVector(n)
# print(a0)

# create using array
a1 = UnifiedVector(n * [1])
# print(a1)

# create using BaseVector
ang = BaseVector(n * [1])
a2 = UnifiedVector(ang)
# print(a1)

# create using UnifiedVector
# the host data is only copied automatically, when the device is not up-to-date!
a2 = UnifiedVector(a1)
a2.UpdateHost()
# print(a2)

# access values
# only changes data on host!
a0[2] = 3
print(a0)

# inner product unified and base
# calculations on host to avoid unnecessary host2device cpy
res = a0.InnerProduct(ang)
print("<a,b>", res)

# inner product unified
res = a0.InnerProduct(a1)
print("innerproduct unified, base:", res)

# scale
a1.Scale(3)
a1.UpdateHost()
# print(a1)

# add to vector
a0.Add(a1, 1)
a0.UpdateHost()
