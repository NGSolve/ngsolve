from ngsolve import BaseMatrix, BitArray, BilinearForm, BaseVector

class SuperLU(BaseMatrix):
    # def __init__(self, a: BaseMatrix | BilinearForm, freedofs: BitArray = None):
    def __init__(self, a, freedofs: BitArray = None):        
        super().__init__()
        self.a = a
        self.freedofs = freedofs

    def Update(self):
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        a = self.a if isinstance(self.a, BaseMatrix) else self.a.mat
        mat = sp.csr_matrix(a.CSR())
        if self.freedofs is not None:
            self.fd = list(self.freedofs)
            mat = mat[self.fd,:][:,self.fd]
        self.lu = spla.factorized(sp.csc_matrix(mat))

    def Mult(self, x: BaseVector, y: BaseVector):
        if not hasattr(self, "lu"):
            self.Update()
        if self.freedofs is not None:
            y.FV().NumPy()[self.fd] = self.lu(x.FV().NumPy()[self.fd])
        else:
            y.FV().NumPy()[:] = self.lu(x.FV().NumPy())
