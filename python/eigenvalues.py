from ngsolve.la import InnerProduct, MultiVector
from math import sqrt
from ngsolve import Projector, Norm, Matrix, Vector, IdentityMatrix

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
    import scipy.linalg

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
        r.SetRandom()
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




def PINVIT(mata, matm, pre, num=1, maxit=20, printrates=True, GramSchmidt=True):
    """preconditioned inverse iteration"""
    import scipy.linalg

    r = mata.CreateRowVector()
    
    uvecs = MultiVector(r, num)
    vecs = MultiVector(r, 2*num)
    # hv = MultiVector(r, 2*num)

    for v in vecs[0:num]:
        v.SetRandom()
    uvecs[:] = pre * vecs[0:num]
    lams = Vector(num * [1])
    
    for i in range(maxit):
        vecs[0:num] = mata * uvecs - (matm * uvecs).Scale (lams)
        vecs[num:2*num] = pre * vecs[0:num]
        vecs[0:num] = uvecs

        vecs.Orthogonalize(matm)

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

        uvecs[:] = vecs * Matrix(evec[:,0:num])
    return lams, uvecs


def LOBPCG(mata, matm, pre, num=1, maxit=20, initial=None, printrates=True, largest=False):
    """Knyazev's cg-like extension of PINVIT"""
    import scipy.linalg

    r = mata.CreateRowVector()

    if initial:
        num=len(initial)
        uvecs = initial
    else:
        uvecs = MultiVector(r, num)

    vecs = MultiVector(r, 3*num)
    for v in vecs:
        r.SetRandom()
        v.data = pre * r

    if initial:
         vecs[0:num] = uvecs       
        
    lams = Vector(num * [1])
    
    for i in range(maxit):
        uvecs.data = mata * vecs[0:num] - (matm * vecs[0:num]).Scale (lams)
        vecs[2*num:3*num] = pre * uvecs
        
        vecs.Orthogonalize(matm)

        asmall = InnerProduct (vecs, mata * vecs)
        msmall = InnerProduct (vecs, matm * vecs)
    
        ev,evec = scipy.linalg.eigh(a=asmall, b=msmall)

        if not largest:
            lams = Vector(ev[0:num])
            if printrates:
                print (i, ":", list(lams), flush=True)

            uvecs[:] = vecs * Matrix(evec[:,0:num])
            vecs[num:2*num] = vecs[0:num]
            vecs[0:num] = uvecs
        else:
            lams = Vector(ev[2*num:3*num])
            if printrates:
                print (i, ":", list(lams), flush=True)

            uvecs[:] = vecs * Matrix(evec[:,2*num:3*num])
            vecs[num:2*num] = vecs[0:num]
            vecs[0:num] = uvecs
        
    return lams, uvecs





def Arnoldi (mat, tol=1e-10, maxiter=200):
    import scipy.linalg
    H = Matrix(maxiter,maxiter, complex=mat.is_complex)
    H[:,:] = 0
    v = mat.CreateVector(colvector=False)
    abv = MultiVector(v, 0)
    v.SetRandom()
    v /= Norm(v)

    for i in range(maxiter):
        abv.Append(v)
        v = (mat*v).Evaluate()
        for j in range(i+1):
            H[j,i] = InnerProduct(v, abv[j])
            v -= H[j,i]*abv[j]
        if i+1 < maxiter:
            H[i+1,i] = Norm(v)
        v = 1/Norm(v)*v

    lam,ev = scipy.linalg.eig(H)
    return Vector(lam), (abv*Matrix(ev)).Evaluate()
    
    





# SOAR: A SECOND-ORDER ARNOLDI METHOD FOR THE SOLUTION OF THE QUADRATIC EIGENVALUE PROBLEM
# Z. Bai and Y. Su, SIAM J. Matrix Anal. Appl 26, pp 640-659 (2005)
# author: A. Schoefl

def SOAR (A, B, maxiter=200):
# first version without memory saving and breakdown 

    q = A.CreateVector()
    p = A.CreateVector()
    s = A.CreateVector()
    r = A.CreateVector()

    Q = MultiVector(q,0)
    P = MultiVector(p,0)
    p[:] = 0
    q.SetRandom()
    q.FV().imag = 0      
    q /= Norm(q)
    T = Matrix(maxiter, complex=A.is_complex)
    T[:,:] = 0
    
    for j in range(maxiter):
        Q.Append(q)
        P.Append(p)
        r.data = A*q + B*p
        s.data = q

        for i in range(j+1):

            T[i,j] = InnerProduct(r, Q[i])
            r -= T[i,j]*Q[i]
            s -= T[i,j]*P[i]
        
        if j+1 < maxiter:
            T[j+1,j] = Norm(r)
            if T[j+1,j] == 0:
                print("SOAR stopped at iteration j = ", j)
                break

            q.data = 1/T[j+1,j]*r
            p.data = 1/T[j+1,j]*s

    return Q




# author: A. Schoefl
def TOAR (A, B, maxiter=200):

    r = A.CreateVector()

    tmp_np = np.zeros((len(r), 2))#, dtype=complex)
    r.SetRandom()
    r.FV().imag = 0      
    r /= Norm(r)
    tmp_np[:,0] = r.FV().real
    r.SetRandom()
    r.FV().imag = 0      
    r /= Norm(r)
    tmp_np[:,1] = r.FV().real

    # tmp_np[:,1] = tmp_np[:,0] # just for testing


    Q_np, X_np, perm = la.qr(tmp_np, pivoting=True, mode="economic")
    # print(Q_np)
    # print(Q_np.shape, X_np.shape, X_np)

    Q = MultiVector(r,0)
    # r.FV().NumPy()[:] = Q_np[:,0]
    
    r.FV().real[:] = Q_np[:,0]
    r.FV().imag = 0
    Q.Append(r)
    # assign rank
    if np.isclose(X_np[1,1], 0):
        eta = 1
    else:
        eta = 2

        # r.FV().NumPy()[:] = Q_np[:,1]
        r.FV().real[:] = Q_np[:,1]
        r.FV().imag = 0
        Q.Append(r)

    # print (Q_np)
    # print (Q[0], Q[1])
    print (InnerProduct(Q,Q))
        
    gamma = np.linalg.norm(tmp_np, ord='fro')
    
    U1 = Matrix(eta,1, True) 
    U1.NumPy()[:,0] = X_np[:eta,1]/gamma

    U2 = Matrix(eta,1, True) 
    U2.NumPy()[:,0] = X_np[:eta,0]/gamma
    print ("X_nb = ", X_np[:eta,:])
    print ("U1 = ", U1)
    print ("U2 = ", U2)

    H = Matrix(maxiter, maxiter-1, True)

    # TODO: would be more efficient in C++
    for j in range(maxiter-1):

        # print(U1, "\n", U2)

        r.data = A*(Q*U1[:,j])+B*(Q*U2[:,j])
        s = Vector(eta, True)
        # MGS: orthogonalize r against Q
        for i in range(eta):
            s[i] = InnerProduct(r, Q[i], conjugate=True) # TODO: order of Q[i], r correct? 
            r-=s[i]*Q[i]
        alpha = InnerProduct(r,r, conjugate=True).real
        
        # MGS 
        for i in range(j):
            # TODO: I have to re-read in proof if should be conjugation makes sense
            # (currently everything is real valued)
            H[i,j] = InnerProduct(s, U1[:,i], conjugate=False) + \
                        InnerProduct( U1[:,j], U2[:,i], conjugate=True)
            s -= H[i,j]*U1[:,i]
            U1[:,j] -= H[i,j]*U2[:,i]

        H[j+1,j] = sqrt(alpha+InnerProduct(s,s,conjugate=True).real+
                    InnerProduct(U1[:,j],U1[:,j], conjugate=True).real)

        alpha = sqrt(alpha)

        # breakdown
        if H[j+1,j] == 0: 
            print("breakdown in iteration ", j)
            return Q

        # deflation
        if alpha == 0:
            print("deflation in iteration ", j)
            tmp = U1
            U1 = Matrix(eta, j+2, True) 
            U1[:,:j+1] = tmp
            U1[:,j+1] = 1/H[j+1,j]*s

            tmp = U2
            U2 = Matrix(eta, j+2, True) 
            U2[:,:j+1] = tmp
            U2[:,j+1] = 1/H[j+1,j]*U1[:,j]
        else:
            # update rank
            eta+=1

            tmp = U2
            U2 = Matrix(eta, j+2, True) 
            U2[:-1,:j+1] = tmp
            U2[:-1,j+1] = 1/H[j+1,j]*U1[:,j]
            U2[eta-1, :] = 0

            tmp = U1
            U1 = Matrix(eta, j+2, True) 
            U1[:-1,:j+1] = tmp[:,:]
            U1[:-1,j+1] = 1/H[j+1,j]*s
            U1[eta-1, j+1] = alpha/H[j+1,j]

            Q.Append(1/alpha*r)                        


    return Q

