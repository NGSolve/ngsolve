# saxpy on the GPU through the ngs_gpu interface
# same script as ngsmetal/examples/gpudevice.py, only the kernel source differs

import numpy as np
from ngsolve.ngstd import *
from ngsolve.ngscuda import *       # registers the cuda backend

dev = GetGPUDevice()
print(dev.name, "- unified memory:", dev.unified_memory,
      "- max threads/group:", dev.max_threads_per_group)

# entry points must be extern "C", otherwise the name is mangled
src = """
extern "C" __global__ void saxpy(const float* x, float* y, float a, int n)
{
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    if (id < n) y[id] = a*x[id] + y[id];
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
