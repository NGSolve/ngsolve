# file: unifiedvector.py
# date: 20.09.2022
#
# testing basic functions for unifiedvectors working on device

from ngsolve.la import *
from ngsolve.ngscuda import *

n = 10

# create vectors using size
print("creating from int")
a0 = UnifiedVector(n)
print(len(a0))
print(a0)

# create using array
print("creating from list")
a1 = UnifiedVector(n * [1])
print(a1)

# create using BaseVector
print("creating from BaseVector")
ang = BaseVector(n * [2])
a2 = UnifiedVector(ang)
print(a1)

# create using UnifiedVector
# the host data is only copied automatically, when the device is not up-to-date!
print("creating from UnifiedVector")
a2 = UnifiedVector(a1)
a2.UpdateHost()
print(a2)

# creating using CreateVector
print("creating using CreateVector")
a3 = a0.CreateVector()
print(a3)

# access values
# only changes data on host!
for i in range(len(a0)):
    a0[i] = 3
print(a0)

print("setting to 2")
a1.FV()[:] = 2
print(a1)

print("setting 4 times 2")
a0.Assign(a1, 4)
print(a0)

# inner product unified and base
# calculations on host to avoid unnecessary host2device cpy
res = a0.InnerProduct(ang)
print("<a,b>", res)

# inner product unified
res = a0.InnerProduct(a1)
print("innerproduct unified, base:", res)

# scale
# a1.Scale(3)
# a1.UpdateHost()
# print(a1)

# add to vector
print("add unified vectors")
a0.Add(a1, 1)
a0.UpdateHost()
print(a0)

print("add unified vectors")
print(a0)
print(a1)
a2 = a0 + a1
print(a2.Evaluate())

print("scale")
a3 = 5 * a0
print(a3.Evaluate())

print("set base using unified")
ang[:] = a0
print(ang)
