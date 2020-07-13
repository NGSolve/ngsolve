from ngsolve.la import InnerProduct, MultiVector
from math import sqrt
from ngsolve import Projector, Norm, Matrix, Vector, IdentityMatrix

try:
    import scipy.linalg
    from scipy import random
except:
    pass


def Orthogonalize (vecs, mat):
    mv = []
    for i in range(len(vecs)):
        for j in range(i):
            vecs[i] -= InnerProduct(vecs[i], mv[j]) * vecs[j]
            
        hv = mat.CreateRowVector()
        hv.data = mat * vecs[i]
        norm = sqrt(InnerProduct(vecs[i], hv))
        vecs[i] *= 1/norm
        hv *= 1/norm
        mv.append (hv)


def PINVIT1(mata, matm, pre, num=1, maxit=20, printrates=True, GramSchmidt=False):
    """preconditioned inverse iteration"""

    r = mata.CreateRowVector()
    Av = mata.CreateRowVector()
    Mv = mata.CreateRowVector()

    uvecs = []
    for i in range(num):
        uvecs.append (mata.CreateRowVector())
    
    vecs = []
    for i in range(2*num):
        vecs.append (mata.CreateRowVector())

    for v in uvecs:
        r.FV().NumPy()[:] = random.rand(len(r.FV()))
        v.data = pre * r

    asmall = Matrix(2*num, 2*num)
    msmall = Matrix(2*num, 2*num)
    lams = num * [1]

    for i in range(maxit):
        
        for j in range(num):
            vecs[j].data = uvecs[j]
            r.data = mata * vecs[j] - lams[j] * matm * vecs[j]
            vecs[num+j].data = pre * r

        if GramSchmidt:
            Orthogonalize (vecs, matm)

        for j in range(2*num):
            Av.data = mata * vecs[j]
            Mv.data = matm * vecs[j]
            for k in range(2*num):
                asmall[j,k] = InnerProduct(Av, vecs[k])
                msmall[j,k] = InnerProduct(Mv, vecs[k])

        ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
        lams[:] = ev[0:num]
        if printrates:
            print (i, ":", lams)
    
        for j in range(num):
            uvecs[j][:] = 0.0
            for k in range(2*num):
                uvecs[j].data += float(evec[k,j]) * vecs[k]

    return lams, uvecs




def PINVIT(mata, matm, pre, num=1, maxit=20, printrates=True, GramSchmidt=False):
    """preconditioned inverse iteration"""

    r = mata.CreateRowVector()
    
    uvecs = MultiVector(r, num)
    vecs = MultiVector(r, 2*num)
    hv = MultiVector(r, 2*num)

    for v in vecs[0:num]:
        v.SetRandom()
    uvecs[:] = pre * vecs[0:num]
    lams = Vector(num * [1])
    
    for i in range(maxit):
        # vecs[0:num] = mata * uvecs
        # vecs[0:num] += (matm * uvecs).Scale (-lams[0:num])

        vecs[0:num] = mata * uvecs - (matm * uvecs).Scale (lams)
        
        # vecs[num:2*num] += (matm * uvecs).Scale (-lams[0:num])
        # for j in range(num):
        # vecs[j] += vecs[num+j]
            # vecs[j] -= lams[j] * vecs[num+j]
            
        vecs[num:2*num] = pre * vecs[0:num]
        vecs[0:num] = uvecs
        
        # vecs[0:num] = uvecs
        # for j in range(num):
        # uvecs[j] = mata * vecs[j] - lams[j] * matm * vecs[j]
        # vecs[num:2*num] = pre * uvecs

        # wish:
        # uvecs[:] = mata * vecs[0:num]  - matm * vecs[0:num] * DiagonalMatrix(lam)

        vecs.Orthogonalize() # matm)

        # hv[:] = mata * vecs
        # asmall = InnerProduct (vecs, hv)
        # hv[:] = matm * vecs
        # msmall = InnerProduct (vecs, hv)
        asmall = InnerProduct (vecs, mata * vecs)
        msmall = InnerProduct (vecs, matm * vecs)
    
        ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)
        lams = Vector(ev[0:num])
        if printrates:
            print (i, ":", list(lams))

            # for j in range(num):
            # uvecs[j] = vecs * Vector(evec[:,j])
        uvecs[:] = vecs * Matrix(evec[:,0:num])
    return lams, uvecs


