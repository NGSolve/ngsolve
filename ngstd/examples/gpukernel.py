# one kernel source, run on the host reference backend and on the gpu,
# and compared - this is the intended shape of a gpu test

import numpy as np
from ngsolve.gpu import *

src = """
KERNEL(saxpy, GLOBAL_IN(float,x), GLOBAL(float,y), VALUE(float,a), VALUE(int,n))
{
  int i = GLOBAL_ID_X;
  if (i < n) y[i] = a*x[i] + y[i];
}

KERNEL(blocksum, GLOBAL_IN(float,x), GLOBAL(float,out), VALUE(int,n))
{
  SHARED(float, s, 256);
  int lid = LOCAL_ID_X;
  int gid = GLOBAL_ID_X;
  s[lid] = (gid < n) ? x[gid] : 0.0f;
  BARRIER();
  for (int stride = 128; stride > 0; stride >>= 1)
    { if (lid < stride) s[lid] += s[lid+stride]; BARRIER(); }
  if (lid == 0) out[GROUP_ID_X] = s[0];
}
"""

n, tg = 1000, 256
ngrp = (n+tg-1)//tg
hx = np.arange(n, dtype=np.float32)

def run(dev):
    lib = dev.CompileSource(GPUKernelPrelude + src)
    q = dev.DefaultQueue()
    x, y, o = dev.NewBuffer(n), dev.NewBuffer(n), dev.NewBuffer(ngrp)
    x.H2D(hx)
    y.H2D(np.ones(n, dtype=np.float32))
    q.Launch(lib.GetKernel("saxpy"), groups=[ngrp], groupsize=[tg], args=[x, y, 2.0, n])
    q.Launch(lib.GetKernel("blocksum"), groups=[ngrp], groupsize=[tg], args=[x, o, n])
    q.Finish()
    return y.D2H(n), o.D2H(ngrp)

cpu, gpu = GetCPUDevice(), GetGPUDevice()
cy, co = run(cpu)
gy, go = run(gpu)

print(f"{cpu.name:15s} saxpy {cy[:4]}  blocksum {co}")
print(f"{gpu.name:15s} saxpy {gy[:4]}  blocksum {go}")
print("saxpy    max diff:", np.abs(cy-gy).max())
print("blocksum max diff:", np.abs(co-go).max())
print("numpy reference  :", (2*hx+1)[:4])
