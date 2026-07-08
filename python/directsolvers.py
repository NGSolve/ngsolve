from ngsolve import BaseMatrix, BitArray, BilinearForm, BaseVector
import ngsolve

class SuperLU(ngsolve.la.SparseFactorizationInterface):
    def Factor(self):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        vals, rows, cols = self.GetInnerMatrix().CSR()
        self.inv_mat = spla.splu(sp.csr_matrix((vals, rows, cols)).tocsc())

    def Solve(self, rhs, sol):
        sol.FV().NumPy()[:] = self.inv_mat.solve(rhs.FV().NumPy())

    def SolveTrans(self, rhs, sol):
        sol.FV().NumPy()[:] = self.inv_mat.solve(rhs.FV().NumPy(), trans="T")

    def SolveConjTrans(self, rhs, sol):
        sol.FV().NumPy()[:] = self.inv_mat.solve(rhs.FV().NumPy(), trans="H")

ngsolve.la.RegisterInverseType("superlu", SuperLU)
