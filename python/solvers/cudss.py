
try:
    import nvmath
    import nvmath.sparse.advanced as nvs
except ImportError:
    raise ImportError("CUDSS solver requires nvmath-python module.")
import ngsolve.la as ngla
import scipy.sparse as sp
import numpy as np

from ngsolve import TimeFunction


class CudssSolver(ngla.SparseFactorizationInterface):
    @TimeFunction
    def Analyze(self):
        self._mat = self.GetInnerMatrix()
        csr = sp.csr_matrix(self._mat.CSR())
        self._tmp = np.empty(csr.shape[1], dtype=csr.dtype)
        options = make_directsolver_options()
        self.solver = nvs.DirectSolver(csr, self._tmp, options=options)
        self.solver.plan()

    @TimeFunction
    def Factor(self):
        self.solver.factorize()

    @TimeFunction
    def Solve(self, b, sol):
        self.solver.reset_operands(b=b.FV().NumPy())
        sol.FV().NumPy()[:] = self.solver.solve()


ngla.RegisterInverseType("cudss", CudssSolver)


# find cudss multithreading lib from installed distribution

import os, pathlib
from importlib import util as importlib_util
from importlib import metadata as importlib_metadata

def _from_dist_files():
    # Use the wheel’s file manifest (most reliable)
    candidates = []
    for dist_name in ("nvidia-cudss-cu12", "nvidia_cudss_cu12", "nvidia-cudss"):  # try common variants
        try:
            dist = importlib_metadata.distribution(dist_name)
        except importlib_metadata.PackageNotFoundError:
            continue
        for f in dist.files or []:
            name = f.name.lower()
            # check for name ends in .so or .so.*
            endwith_so = name.endswith(".so") or (".so." in name and name.rsplit(".so.", 1)[1].replace(".", "").isdigit())
            if name.startswith("libcudss_mtlayer_") and endwith_so:
                candidates.append(dist.locate_file(f))
            if name.startswith("cudss_mtlayer_") and name.endswith(".dll"):
                candidates.append(dist.locate_file(f))
        if candidates:
            # Prefer anything in bin/ or lib/ if multiple
            candidates.sort(key=lambda p: ("bin" not in str(p) and "lib" not in str(p), str(p)))
            return str(candidates[0])
    return None


def make_directsolver_options():
    # Helpful on Windows (Python 3.8+): ensure DLL deps can be found
    if os.name == "nt":
        for var in ("CUDA_PATH", "CUDSS_PATH", "CONDA_PREFIX"):
            base = os.environ.get(var)
            if base:
                p = pathlib.Path(base) / "bin"
                if p.exists():
                    try:
                        os.add_dll_directory(str(p))
                    except Exception:
                        pass
    mtlib = _from_dist_files()
    return nvs.DirectSolverOptions(multithreading_lib=mtlib)
