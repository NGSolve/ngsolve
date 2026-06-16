import numpy as np
import ctypes
import ctypes.util

from .. import la as ngla
from .. import config, TimeFunction

_mkl_free_buffers = None
_pardiso = None


def _find_mkl():
    # Register the solvers to provide a proper error message if MKL Pardiso is not available
    ngla.RegisterInverseType("pardiso", MKLPardiso)
    ngla.RegisterInverseType("pardisospd", MKLPardisoSPD)

    import importlib.metadata

    global _mkl_free_buffers, _pardiso

    if _mkl_free_buffers is not None and _pardiso is not None:
        return

    mkl_from_pip = False
    libname = None
    try:
        importlib.metadata.version("mkl")
    except importlib.metadata.PackageNotFoundError:
        return

    # iterate over all files in the mkl package
    dist = importlib.metadata.distribution("mkl")
    for file in dist.files or []:
        filename = str(dist.locate_file(file))
        if "mkl_rt" in filename:
            libname = filename
            mkl_from_pip = True
            break
    if libname is None:
        libname = ctypes.util.find_library("mkl_rt")
    if libname is None:
        return
    try:
        mkl_rt = ctypes.CDLL(libname)
    except OSError:
        return

    try:
        _mkl_free_buffers = mkl_rt.mkl_free_buffers
        _mkl_free_buffers.restype = None
        _mkl_free_buffers.argtypes = []
    except AttributeError:
        _mkl_free_buffers = None

    try:
        _pardiso = mkl_rt.pardiso
        _pardiso.restype = None
    except AttributeError:
        _pardiso = None

    if _pardiso is not None and _mkl_free_buffers is not None and mkl_from_pip:
        ngla.BaseMatrix.SetDefaultInverseType("pardiso")


class MKLPardiso(ngla.SparseFactorizationInterface):
    is_spd = False

    @staticmethod
    def _as_ctypes(arr):
        if arr.dtype == np.complex128:
            return np.ctypeslib.as_ctypes(arr.view(np.float64))
        return np.ctypeslib.as_ctypes(arr)

    def __init__(
        self, *args, msglevel=0, params: dict[int, int] | None = None, **kwargs
    ):
        from ctypes import c_int

        if _pardiso is None or _mkl_free_buffers is None:
            raise RuntimeError(
                "MKL Pardiso is not available. Ensure that MKL is installed: pip install mkl"
            )

        super().__init__(*args, **kwargs)

        self._params = (c_int * 64)()
        for i in range(64):
            self._params[i] = 0

        if params is not None:
            for key, value in params.items():
                if 0 <= key < 64:
                    self._params[key] = value

        self._params[0] = 1  # no pardiso defaults
        self._params[1] = 0  # fill in 0..MDO, 2..metis
        self._params[2] = 8  # nthreads
        self._params[7] = 16
        self._params[9] = 13  # perturbation 1E-10
        self._params[10] = 1
        self._params[12] = 1  # slicing + matching
        self._params[17] = -1
        self._params[20] = 1  # 1x1 and 2x2 Bunch-Kaufman pivoting
        self._params[26] = 1  # check input matrix
        self._params[34] = 1  # 0-based (C-style) indexing
        self._params[59] = 0  # 0..incore

        self._pt = (c_int * 128)()
        for i in range(128):
            self._pt[i] = 0

        int_arg_t = c_int * 1

        self._maxfct = int_arg_t(1)
        self._mnum = int_arg_t(1)
        self._matrixtype = int_arg_t(11)
        self._phase = int_arg_t(12)
        self._nrhs = int_arg_t(1)
        self._msglevel = int_arg_t(msglevel)
        self._error = int_arg_t(0)
        self._n = int_arg_t(0)

    def _call_pardiso(self, phase):
        self._phase[0] = phase
        null = ctypes.c_void_p(0)
        _pardiso(
            self._pt,
            self._maxfct,
            self._mnum,
            self._matrixtype,
            self._phase,
            self._n,
            self._data,
            self._indptr,
            self._indices,
            null,
            self._nrhs,
            self._params,
            self._msglevel,
            self._b,
            self._x,
            self._error,
        )
        if self._error[0] != 0:
            raise RuntimeError(f"MKL Pardiso error in phase {phase}: {self._error[0]}")

    def get_csr(self):
        mat = self.GetInnerMatrix()
        if self.is_symmetric.is_true:
            mat = (
                mat.CreateTranspose()
                if self.is_symmetric_storage
                else ngla.ExtractTri(mat, lower=False)
            )
            data = [np.array(v.NumPy()) for v in mat.CSR()]
        else:
            data = [v.NumPy() for v in mat.CSR()]
        data[2] = np.asarray(data[2], dtype=np.int32)
        return data, mat

    @TimeFunction
    def Analyze(self):
        if hasattr(self, "_pt"):
            self._release()

        data, self._csr_mat = self.get_csr()
        self._data_np = np.array(data[0])  # writable copy for pardiso
        self._indices_np = data[1]  # view kept alive by self._csr_mat
        self._indptr_np = data[2]  # already a copy (dtype conversion)

        is_symmetric = self.is_symmetric.is_true
        self._params[12] = 0 if is_symmetric else 1  # slicing + matching

        if self.is_complex:
            self._matrixtype[0] = 6 if is_symmetric else 13
        else:
            if self.is_spd:
                self._matrixtype[0] = 2
            else:
                self._matrixtype[0] = -2 if is_symmetric else 11

        self._n[0] = n = self.GetInnerMatrix().height

        self._data = self._as_ctypes(self._data_np)
        self._indices = self._as_ctypes(self._indices_np)
        self._indptr = self._as_ctypes(self._indptr_np)

        dtype = np.complex128 if self.is_complex else np.float64
        self._b_np = np.zeros(n, dtype=dtype)
        self._x_np = np.zeros(n, dtype=dtype)
        self._b = self._as_ctypes(self._b_np)
        self._x = self._as_ctypes(self._x_np)

        self._call_pardiso(12)
        self._is_first_factor_call = True

    @TimeFunction
    def Factor(self):
        if not self._is_first_factor_call:
            data, _ = self.get_csr()
            self._data_np[:] = data[0]
            self._call_pardiso(22)
        self._is_first_factor_call = False

    @TimeFunction
    def Solve(self, b, sol):
        self._b_np[:] = b.FV().NumPy()

        self._call_pardiso(33)

        sol.FV().NumPy()[:] = self._x_np

    def _release(self):
        if hasattr(self, "_pt"):
            try:
                self._phase[0] = -1
                null = ctypes.c_void_p(0)
                _pardiso(
                    self._pt,
                    self._maxfct,
                    self._mnum,
                    self._matrixtype,
                    self._phase,
                    self._n,
                    null,
                    self._indptr,
                    self._indices,
                    null,
                    self._nrhs,
                    self._params,
                    self._msglevel,
                    null,
                    null,
                    self._error,
                )
            except Exception:
                pass

    def __del__(self):
        self._release()
        if _mkl_free_buffers is not None:
            _mkl_free_buffers()


class MKLPardisoSPD(MKLPardiso):
    is_spd = True


if not config.USE_MKL and not config.USE_PARDISO:
    # no pardiso was linked to NGSolve, try to find it at run time
    _find_mkl()
