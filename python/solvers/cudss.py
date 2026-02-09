try:
    import nvmath.sparse.advanced as nvs
    from nvmath.internal import utils
    from nvmath.internal.tensor_ifc_numpy import NumpyTensor
except ImportError:
    raise ImportError("CUDSS solver requires nvmath-python module.")

import ngsolve.la as ngla
import scipy.sparse as sp
import numpy as np

from ngsolve import TimeFunction

class CudssSolver(ngla.SparseFactorizationInterface):
    solver: nvs.DirectSolver = None

    @TimeFunction
    def Analyze(self):
        if hasattr(self, "solver") and self.solver is not None:
            self.solver.free()
            del self.solver

        csr = sp.csr_matrix(self.GetInnerMatrix().CSR())

        options = make_directsolver_options()
        self.extract_symmetric = self.is_symmetric.is_true or (csr != csr.T).nnz == 0
        if self.extract_symmetric:
            if not self.is_symmetric_storage:
                csr = sp.tril(csr, format="csr")
            options.sparse_system_type = nvs.DirectSolverMatrixType.SYMMETRIC
            options.sparse_system_view = nvs.DirectSolverMatrixViewType.LOWER

        tmp = np.empty(csr.shape[1], dtype=csr.dtype)
        self.solver = nvs.DirectSolver(csr, tmp, options=options)
        self.solver.plan()
        self._is_first_factor_call = True

    @TimeFunction
    def Factor(self):
        if not self._is_first_factor_call:
            mat = self.GetInnerMatrix()
            if self.extract_symmetric and not self.is_symmetric_storage:
                values = sp.tril(sp.csr_matrix(mat.CSR()), format="csr").data
            else:
                values = mat.AsVector().FV().NumPy()
            stream_holder = utils.get_or_create_stream(self.solver.device_id, None, self.solver.rhs_package)
            values_tensor = NumpyTensor(values)
            self.solver.a.values.copy_(values_tensor, stream_holder)
        self.solver.factorize()
        self._is_first_factor_call = False

    @TimeFunction
    def Solve(self, b, sol):
        stream_holder = utils.get_or_create_stream(self.solver.device_id, None, self.solver.rhs_package)
        self.solver.b.copy_(NumpyTensor(b.FV().NumPy()), stream_holder)
        sol.FV().NumPy()[:] = self.solver.solve()

    def __del__(self):
        if self.solver is not None:
            self.solver.free()
            del self.solver


ngla.RegisterInverseType("cudss", CudssSolver)


# find cudss multithreading lib from installed distribution

import os, pathlib
from importlib import util as importlib_util
from importlib import metadata as importlib_metadata

def _from_dist_files():
    # Use the wheel’s file manifest (most reliable)
    candidates = []
    # for dist_name in ("nvidia-cudss-cu12", "nvidia_cudss_cu12", "nvidia-cudss"):  # try common variants
    for dist_name in ("nvidia-cudss-cu12", "nvidia_cudss_cu12", "nvidia-cudss",
                      "nvidia-cudss-cu13", "nvidia_cudss_cu13"):  # try common variants
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
