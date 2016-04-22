from ngsolve.la import InnerProduct
from math import sqrt

def pcg(mat, pre, rhs, maxits = 100):
    """preconditioned conjugate gradient method"""

    u = rhs.CreateVector()
    d = rhs.CreateVector()
    w = rhs.CreateVector()
    s = rhs.CreateVector()

    u[:] = 0.0
    d.data = rhs - mat * u
    w.data = pre * d
    s.data = w
    wdn = InnerProduct (w,d)
    print ("|u0| = ", wdn)
    
    for it in range(maxits):
        w.data = mat * s
        wd = wdn
        Ass = InnerProduct (s, w)
        alpha = wd / Ass
        u.data += alpha * s
        d.data += (-alpha) * w

        w.data = pre * d
        wdn = InnerProduct (w, d)
        beta = wdn / wd

        s *= beta
        s.data += w

        err = sqrt(wd)
        if err < 1e-16:
            break
        print ("it = ", it, " err = ", err)

    return u




