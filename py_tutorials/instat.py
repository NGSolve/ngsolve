from ngsolve.comp import PyNumProc
from ngsolve.solve import Redraw

from math import sin
from time import sleep


class npParabolic(PyNumProc):

    def __init__(self, pde, flags):
        super(npParabolic,self).__init__(pde,flags)
        self.mypde = pde

        
    def Do(self, heap):
        print ("solve parabolic equation")
            
        pde = self.pde
        
        tau = 0.1
        
        v = pde.spaces["v"]
        u = pde.gridfunctions["u"].vec
        f = pde.linearforms["f"].vec
        a = pde.bilinearforms["a"].mat
        m = pde.bilinearforms["m"].mat
        
        hmat = a.CreateMatrix()
        hmat.AsVector().data = tau * a.AsVector() + m.AsVector()
        
        inv = hmat.Inverse(v.FreeDofs())
        
        d = u.CreateVector()
        w = u.CreateVector()
        
        t = 0.0;
        for j in range (0,100000):
            t += tau
            d.data = sin(t) * f - a * u
            w.data = inv * d
            u.data += tau * w
            
            print ("t = ", t)
            Redraw()
            sleep (0.05)
           
       
   
