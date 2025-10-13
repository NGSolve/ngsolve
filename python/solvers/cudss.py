
try:
    import nvmath
    import nvmath.sparse.advanced as nvs
except ImportError:
    raise ImportError("CUDSS solver requires nvmath-python module.")
import ngsolve.la as ngla
import scipy.sparse as sp
import numpy as np


class CudssSolver(ngla.SparseFactorizationInterface):
    def Analyze(self):
        self._mat = self.GetInnerMatrix()
        csr = sp.csr_matrix(self._mat.CSR())
        self._tmp = np.empty(csr.shape[1], dtype=csr.dtype)
        self.solver = nvs.DirectSolver(csr, self._tmp)
        self.solver.plan()

    def Factor(self):
        self.solver.factorize()

    def Solve(self, b, sol):
        self.solver.reset_operands(b=b.FV().NumPy())
        sol.FV().NumPy()[:] = self.solver.solve()


ngla.RegisterInverseType("cudss", CudssSolver)
