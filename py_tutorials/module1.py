print ("**** this is my first ngs-py module *******")

from ngsolve import solve

def Hi():
     print ("hello")


class pyNP1(solve.PyNumProc):

   def __init__(self, pde, flags):
        print ('********* pyNP1 constructor ********')
        super(pyNP1,self).__init__(pde,flags)

   def Do(self, heap):
        print('******** pyNP1::Do ********************')
        print (self.pde.gridfunctions["u"])
        u = self.pde.gridfunctions["u"]
        print ("ndof = ", u.space.ndof)
        print('*************************************')



print ("*********** module done ****************")



