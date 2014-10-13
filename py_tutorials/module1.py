print ("**** this is my first ngsolve-py module *****")

from ngsolve.comp import PyNumProc
from cg import *
from ngsolve.solve import Redraw


class pyNP1(PyNumProc):

   def __init__(self, pde, flags):
        print ('********* pyNP1 constructor ********')
        super(pyNP1,self).__init__(pde,flags)

   def Do(self, heap):
        print('******** pyNP1::Do ********************')
        print (self.pde.gridfunctions["u"])
        u = self.pde.gridfunctions["u"]
        print ("ndof = ", u.space.ndof)
        print('*************************************')




class pybvp(PyNumProc):
   def __init__(self, pde, flags):
        super(pybvp,self).__init__(pde,flags)
        self.mypde = pde
        
   def Do(self, heap):
        print('******** pybvp::Do ********************')

        A = self.mypde.bilinearforms["a"]
        C = self.mypde.preconditioners["c"]
        f = self.mypde.linearforms["f"]
        u = self.mypde.gridfunctions["u"]

        u.vec.data = pcg (mat = A.mat, pre = C.mat, rhs = f.vec, maxits = 50)
        Redraw()
    




print ("*********** module done ****************")



