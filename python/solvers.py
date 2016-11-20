from ngsolve.la import InnerProduct
from math import sqrt


def CG(mat, rhs, pre=None, sol=None, maxsteps = 100, printrates = True, initialize = True):
    """preconditioned conjugate gradient method"""

    u = sol if sol else rhs.CreateVector()
    d = rhs.CreateVector()
    w = rhs.CreateVector()
    s = rhs.CreateVector()

    if initialize: u[:] = 0.0
    d.data = rhs - mat * u
    w.data = pre * d if pre else d
    s.data = w
    wdn = InnerProduct (w,d)
    
    for it in range(maxsteps):
        w.data = mat * s
        wd = wdn
        as_s = InnerProduct (s, w)
        alpha = wd / as_s
        u.data += alpha * s
        d.data += (-alpha) * w

        w.data = pre*d if pre else d
        
        wdn = InnerProduct (w, d)
        beta = wdn / wd

        s *= beta
        s.data += w

        err = sqrt(wd)
        if err < 1e-16: break
            
        if printrates:
            print ("it = ", it, " err = ", err)

    return u



