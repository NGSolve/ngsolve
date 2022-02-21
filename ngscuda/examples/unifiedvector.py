# file: unifiedvector.py
# date: 28.09.2022
#
# testing basic functions for unifiedvectors working on device
# UnifiedVector stores host and device data

from ngsolve.la import *
from ngsolve.ngscuda import *

n = 3

# Create UnifiedVector using size
print("creating from size")
a0 = UnifiedVector(n)
# print(a0)

# Create UnifiedVector using array
print("creating from list")
a1 = UnifiedVector(n * [1])
# print(a1)

# Create UnifiedVector using BaseVector
print("creating from BaseVector")
ang = BaseVector(n)
a2 = UnifiedVector(ang)
# print(a1)

# Create UnifiedVector using UnifiedVector
# the host data is only copied, when the device is not up-to-date!
print("creating from UnifiedVector")
a2 = UnifiedVector(a1)
a2.UpdateHost()
# print(a2)

# Create UnifiedVector using CreateVector
# values are undefined!
print("creating using CreateVector")
a3 = a1.CreateVector()
# print(a3)



# Set values using scalar
print("setting to 2")
a0.FV()[:] = 1
a1.FV()[:] = 2
a2.FV()[:] = 3
# print(a1)

# Set BaseVector using UnifiedVector
print("set basevector using unifiedvector")
ang[:] = a1
# print(ang)



# inner product unified and base
# calculations on host to avoid unnecessary host2device cpy
res = a0.InnerProduct(ang)
# print("innerproduct unified, base:", res)

# inner product unified
res = a0.InnerProduct(a1)
# print("innerproduct unified, unified:", res)

# add unifed and unified
print("add unified vectors")
a2 = a0 + a1
# print(a2.Evaluate())

# scale unified
print("scale")
a3 = 5 * a0
# print(a3.Evaluate())

