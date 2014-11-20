from ngsolve.comp import PyNumProc
from ngsolve.fem import *



class npALE(PyNumProc):
        
    def Do(self, heap):
        print ("ALE - Laplace")

        pde = self.pde

        d = pde.gridfunctions["def"]
        d.Set(VariableCF("(2*y*(1-y)*x*(1-x),(3*y*(1-y)*x*(1-x))"))
        pde.Mesh().SetDeformation (d)

        u = pde.gridfunctions["u"]
        # first test:
        # u.Set(VariableCF("x"))

        
        v = pde.spaces["v"]
        a = pde.bilinearforms["a"]
        f = pde.linearforms["f"]

        a.ReAssemble()
        f.ReAssemble()
        inv = a.mat.Inverse(v.FreeDofs())

        u.vec.data = inv * f.vec
        
