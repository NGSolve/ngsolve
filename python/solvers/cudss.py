import ngsolve.la as ngla

import numpy as np
from ngsolve import TimeFunction

def _stream_package(solver):
    # nvmath >= 1.0 dropped "numpy" from the stream package registry, leaving "cuda".
    from nvmath.internal import package_wrapper
    pkg = solver.rhs_package
    return pkg if pkg in package_wrapper.PACKAGE else "cuda"


def _free_device_memory():
    from cuda.bindings import runtime
    err, free, total = runtime.cudaMemGetInfo()
    return free


def _set_hybrid_limit(solver, limit=None):
    # keep the factors in host memory, limit is the device memory used as working space
    # (bytes, "6GiB", "80%"), default 90% of the device memory left after the analysis
    from nvmath.internal import utils
    free = _free_device_memory()
    need = int(solver.plan_info.memory_estimates.hybrid_min_device_memory)
    if limit is None:
        limit = max(int(0.9 * free), min(need, free))
    else:
        limit = utils.get_memory_limit_from_device_id(int(limit) if isinstance(limit, float) else limit, solver.device_id)
    if limit < need:
        raise MemoryError(f"cuDSS hybrid memory mode needs {need / 2**30:.2f} GiB device memory, "
                          f"the limit is {limit / 2**30:.2f} GiB ({free / 2**30:.2f} GiB free)")
    solver._internal_config.hybrid_device_memory_limit = limit


class CudssSolver(ngla.SparseFactorizationInterface):
    solver = None
    dtype = np.float64
    spd = False

    def __init__(self, *args, **kwargs):
        """inverse flags:
        hybrid: True or device memory limit (bytes, "6GiB", "80%") to keep the factors
          in host memory, False to never do so. Default: switch to hybrid memory mode
          if the factors do not fit into device memory.
        reordering: default|nested_dissection|amd|colamd|btf_colamd. colamd enables
          global partial pivoting (more fill, no hybrid memory mode).
        pivot_zero_diagonal: cuDSS pivots only locally by default and fails on zero
          diagonal entries (saddle point systems). They are perturbed by this factor
          times the row maximum and the solution is corrected by iterative refinement
          with the original matrix. Default 1e-8, False to disable.
        refinement: maximal iterative refinement steps, default 3 if a diagonal is
          perturbed, else 0."""
        super().__init__(*args, **kwargs)
        flags = {}
        if args and hasattr(args[0], "inverse_flags"):
            flags = args[0].inverse_flags.ToDict()
        self.hybrid = flags.get("hybrid")
        self.reordering = flags.get("reordering")
        self.pivot_zero_diagonal = flags.get("pivot_zero_diagonal", 1e-8)
        self.refinement = flags.get("refinement")
        self._hybrid_active = self.hybrid not in (None, False)

    @TimeFunction
    def Analyze(self):
        try:
            import nvmath.sparse.advanced as nvs
        except ImportError:
            raise ImportError("CUDSS solver requires nvmath-python module.")
        try:
            import scipy.sparse as sp
        except ImportError:
            raise ImportError("CUDSS solver requires scipy.")
        if hasattr(self, "solver") and self.solver is not None:
            self.solver.free()
            del self.solver

        inner = self.GetInnerMatrix()
        csr = sp.csr_matrix(inner.CSR())
        if np.iscomplexobj(csr.data):
            self.dtype = np.complex64 if self.dtype == np.float32 else np.complex128
        csr = csr.astype(self.dtype, copy=False)
        self._full = csr

        options = make_directsolver_options()
        use_spd = self.spd or getattr(self, "is_spd", False)
        if use_spd or self.is_symmetric.is_true:
            self.extract_symmetric = True
        else:
            self.extract_symmetric = (csr != csr.T).nnz == 0
        if self.extract_symmetric:
            if not self.is_symmetric_storage:
                csr = sp.tril(csr, format="csr")
            if use_spd:
                options.sparse_system_type = (nvs.DirectSolverMatrixType.HPD
                                              if np.iscomplexobj(csr.data)
                                              else nvs.DirectSolverMatrixType.SPD)
            else:
                options.sparse_system_type = nvs.DirectSolverMatrixType.SYMMETRIC
            options.sparse_system_view = nvs.DirectSolverMatrixViewType.LOWER

        self._diag_pos, self._diag_pert = None, None
        if self.pivot_zero_diagonal and not use_spd:
            zero = abs(self._full.diagonal()) <= 1e-14 * abs(self._full.data).max()
            if zero.any():
                rows = np.repeat(np.arange(csr.shape[0]), np.diff(csr.indptr))
                self._diag_pos = np.where((rows == csr.indices) & zero[rows])[0]
                rowmax = abs(self._full).max(axis=1).toarray().ravel()
                self._diag_pert = self.pivot_zero_diagonal * rowmax[rows[self._diag_pos]]
                csr = csr.copy()
                csr.data[self._diag_pos] += self._diag_pert
        self._nrefine = (self.refinement if self.refinement is not None
                         else 3 if self._diag_pos is not None else 0)

        reordering = self.reordering
        tmp = np.empty(csr.shape[1], dtype=csr.dtype)
        execution = None
        if self._hybrid_active:
            execution = nvs.ExecutionCUDA(hybrid_memory_mode_options=nvs.HybridMemoryModeOptions(hybrid_memory_mode=True))
        self.solver = nvs.DirectSolver(csr, tmp, options=options, execution=execution)
        self._global_pivoting = reordering in ("colamd", "btf_colamd")
        if reordering not in (None, "default"):
            from nvmath.bindings import cudss
            self.solver.plan_config.reordering_algorithm = cudss.ReorderingAlg[reordering.upper()]
        try:
            self.solver.plan()
        except Exception as e:
            from nvmath.bindings.cudss import cuDSSError, Status
            if not (isinstance(e, cuDSSError) and e.status == Status.ALLOC_FAILED):
                raise
            raise MemoryError(f"cuDSS analysis does not fit into device memory ({_free_device_memory() / 2**30:.2f} GiB free)") from e
        self._is_first_factor_call = True
        if self._hybrid_active:
            # the minimal device memory is known after the analysis
            _set_hybrid_limit(self.solver, None if self.hybrid in (None, True) else self.hybrid)
        elif self.hybrid is None and self.solver.plan_info.memory_estimates.peak_device_memory > _free_device_memory():
            self._switch_to_hybrid()

    def _switch_to_hybrid(self):
        if self._global_pivoting:
            raise MemoryError("cuDSS: factors with global pivoting (colamd) do not fit into "
                              "device memory, and hybrid memory mode does not support it")
        print("cuDSS: factors do not fit into device memory, switching to hybrid memory mode")
        self._hybrid_active = True
        self.solver.free()
        self.solver = None
        self.Analyze()

    @TimeFunction
    def Factor(self):
        if not self._is_first_factor_call:
            from nvmath.internal.tensor_ifc_numpy import NumpyTensor
            from nvmath.internal import utils
            import scipy.sparse as sp
            mat = self.GetInnerMatrix()
            if self.extract_symmetric and not self.is_symmetric_storage:
                values = sp.tril(sp.csr_matrix(mat.CSR()), format="csr").data
            else:
                values = mat.AsVector().FV().NumPy()
            values = values.astype(self.dtype, copy=False)
            if self._diag_pos is not None:
                values = values.copy()
                values[self._diag_pos] += self._diag_pert
            stream_holder = utils.get_or_create_stream(self.solver.device_id, None, _stream_package(self.solver))
            values_tensor = NumpyTensor(values)
            self.solver.a.values.copy_(values_tensor, stream_holder)
        try:
            self.solver.factorize()
        except Exception as e:
            from nvmath.bindings.cudss import cuDSSError, Status
            if self._hybrid_active or self.hybrid is not None or not (isinstance(e, cuDSSError) and e.status == Status.ALLOC_FAILED):
                raise
            self._switch_to_hybrid()
            self.solver.factorize()
        self._is_first_factor_call = False

    @TimeFunction
    def Solve(self, b, sol):
        from nvmath.internal import utils
        from nvmath.internal.tensor_ifc_numpy import NumpyTensor
        stream_holder = utils.get_or_create_stream(self.solver.device_id, None, _stream_package(self.solver))
        rhs = b.FV().NumPy().astype(self.dtype, copy=False)
        self.solver.b.copy_(NumpyTensor(rhs), stream_holder)
        x = self.solver.solve().copy()
        if self._nrefine:
            # refinement with the original matrix, the factors may be perturbed
            A = self._full
            bnorm = np.linalg.norm(rhs)
            rnorm_old = np.inf
            for _ in range(self._nrefine):
                r = rhs - A @ x
                rnorm = np.linalg.norm(r)
                if rnorm <= 1e-15 * bnorm or rnorm > 0.5 * rnorm_old:
                    break
                rnorm_old = rnorm
                self.solver.b.copy_(NumpyTensor(r), stream_holder)
                x += self.solver.solve()
        sol.FV().NumPy()[:] = x

    def __del__(self):
        if self.solver is not None:
            self.solver.free()
            del self.solver


class CudssSolverFloat32(CudssSolver):
    dtype = np.float32


class CudssSolverSPD(CudssSolver):
    spd = True


class CudssSolverFloat32SPD(CudssSolver):
    dtype = np.float32
    spd = True


ngla.RegisterInverseType("cudss", CudssSolver)
ngla.RegisterInverseType("cudss_float32", CudssSolverFloat32)
ngla.RegisterInverseType("cudss_spd", CudssSolverSPD)
ngla.RegisterInverseType("cudss_float32_spd", CudssSolverFloat32SPD)


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
    import nvmath.sparse.advanced as nvs
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
