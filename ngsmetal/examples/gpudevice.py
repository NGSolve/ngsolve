# saxpy on the GPU through the ngs_gpu interface

import numpy as np
from ngsolve.ngstd import *
from ngsolve.ngsmetal import *        # registers the metal backend

dev = GetGPUDevice()
print(dev.name, "- unified memory:", dev.unified_memory,
      "- max threads/group:", dev.max_threads_per_group)

src = """
#include <metal_stdlib>
using namespace metal;
kernel void saxpy(device const float* x [[buffer(0)]],
                  device float* y       [[buffer(1)]],
                  constant float& a     [[buffer(2)]],
                  constant int& n       [[buffer(3)]],
                  uint id [[thread_position_in_grid]]) {
    if (id < uint(n)) y[id] = a*x[id] + y[id];
}
"""

n, tg = 1000, 256
kernel = dev.CompileSource(src).GetKernel("saxpy")
queue = dev.DefaultQueue()

x = dev.NewBuffer(n, np.float32)
y = dev.NewBuffer(n, np.float32)
x.H2D(np.arange(n, dtype=np.float32))
y.H2D(np.ones(n, dtype=np.float32))

queue.Launch(kernel, groups=[(n+tg-1)//tg], groupsize=[tg], args=[x, y, 2.0, n])
queue.Finish()

print(y.D2H(n)[:8])        # 2*i + 1
