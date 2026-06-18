import ctypes
import importlib.resources as r

from ngsolve import config as _config

if not _config.USE_CUDA:
    raise ImportError("ngscuda was imported, but CUDA support is not enabled in this build of NGSolve.")

def preload_cuda_libs():
    libs = [
        ("nvidia.cuda_runtime", "lib/libcudart.so.12"),
        ("nvidia.nvjitlink", "lib/libnvJitLink.so.12"),
        ("nvidia.cublas", "lib/libcublas.so.12"),
    ]

    for pkg, rel in libs:
        path = r.files(pkg).joinpath(rel)
        ctypes.CDLL(str(path), mode=ctypes.RTLD_GLOBAL)


if _config.is_python_package:
    preload_cuda_libs()

from ngsolve._ngscuda import *
