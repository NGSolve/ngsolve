"""
Common GPU layer - one import for whichever backend this build has.

    from ngsolve.gpu import *

    dev = GetGPUDevice()
    lib = dev.CompileSource(GPUKernelPrelude + src)

Importing this module registers ngsmetal or ngscuda if one of them is available,
and otherwise selects the host reference device, so GetGPUDevice() never returns
None. Which one was taken is in `backend`: "metal", "cuda" or "host".

Everything the backend module exports (device vectors, device matrices, ...) is
re-exported here as well.
"""

from ngsolve.ngstd import (GPUDevice, GPUBuffer, GPUKernel, GPULibrary, GPUQueue,
                           MemType, GPUKernelPrelude, TinyBlaPrelude,
                           GetGPUDevice, GetCPUDevice, HasGPUDevice, SetGPUDevice)

backend = None
_backend_names = []

for _tag, _module in (("metal", "ngsolve.ngsmetal"), ("cuda", "ngsolve.ngscuda")):
    try:
        _mod = __import__(_module, fromlist=["*"])
    except ImportError:
        continue          # not built, or no runtime - try the next one
    backend = _tag
    for _name in dir(_mod):
        if not _name.startswith("_"):
            globals()[_name] = getattr(_mod, _name)
            _backend_names.append(_name)
    break

# a backend may be built but find no device, so ask rather than assume
if GetGPUDevice() is None:
    SetGPUDevice(GetCPUDevice())
    backend = "host"

__all__ = ["GPUDevice", "GPUBuffer", "GPUKernel", "GPULibrary", "GPUQueue",
           "MemType", "GPUKernelPrelude", "TinyBlaPrelude",
           "GetGPUDevice", "GetCPUDevice", "HasGPUDevice", "SetGPUDevice",
           "backend"] + _backend_names
