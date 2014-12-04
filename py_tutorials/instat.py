from ngsolve.comp import PyNumProc
from ngsolve.solve import Redraw

from math import sin
from time import sleep


class npParabolic(PyNumProc):
        
    def Do(self, heap):
        print ("solve parabolic equation")
        
        tau = 0.1
        
        pde = self.pde
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
            Redraw(blocking=True)
            # sleep (0.001)
           
       
   
