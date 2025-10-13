from ngsolve import BaseMatrix, BitArray, BilinearForm, BaseVector
import ngsolve

class SuperLU(ngsolve.la.SparseFactorizationInterface):
    def Factor(self):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        vals, rows, cols = self.GetInnerMatrix().CSR()
        self.inv_mat = spla.factorized(sp.csr_matrix((vals, rows, cols)))

    def Solve(self, rhs, sol):
        sol.FV().NumPy()[:] = self.inv_mat(rhs.FV().NumPy())

ngsolve.la.RegisterInverseType("superlu", SuperLU)
