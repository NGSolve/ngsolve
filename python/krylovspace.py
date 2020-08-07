
from ngsolve import Projector, Norm, TimeFunction, BaseMatrix, Preconditioner, InnerProduct, \
    Norm, sqrt, Vector, Matrix, BaseVector, BitArray
from typing import Optional, Callable
import logging
from netgen.libngpy._meshing import _PushStatus, _GetStatus, _SetThreadPercentage
from math import log

class CGSolver(BaseMatrix):
    def __init__(self, mat : BaseMatrix, pre : Optional[Preconditioner] = None,
                 freedofs : Optional[BitArray] = None,
                 conjugate : bool = False, tol : float = 1e-12, maxsteps : int = 100,
                 callback : Optional[Callable[[int, float], None]] = None,
                 printing=False, abstol=None):
        super().__init__()
        self.mat = mat
        assert (freedofs is None) != (pre is None) # either pre or freedofs must be given
        self.pre = pre if pre else Projector(freedofs, True)
        self.conjugate = conjugate
        self.tol = tol
        self.abstol = abstol
        self.maxsteps = maxsteps
        self.callback = callback
        self._tmp_vecs = [self.mat.CreateRowVector() for i in range(3)]

        self.printing = printing
        self.errors = []
        self.iterations = 0

    def Height(self) -> int:
        return self.mat.width

    def Width(self) -> int:
        return self.mat.height

    def IsComplex(self) -> bool:
        return self.mat.IsComplex()

    def Mult(self, x : BaseVector, y : BaseVector) -> None:
        self.Solve(rhs=x, sol=y, initialize=True)

    def Update(self):
        if hasattr(self.pre, "Update"):
            self.pre.Update()

    @TimeFunction
    def Solve(self, rhs : BaseVector, sol : Optional[BaseVector] = None,
              initialize : bool = True) -> None:
        old_status = _GetStatus()
        _PushStatus("CG Solve")
        _SetThreadPercentage(0)
        self.sol = sol if sol is not None else self.mat.CreateRowVector()
        d, w, s = self._tmp_vecs
        u, mat, pre, conjugate, tol, maxsteps, callback = self.sol, self.mat, self.pre, self.conjugate, \
            self.tol, self.maxsteps, self.callback
        if initialize:
            u[:] = 0
            d.data = rhs
        else:
            d.data = rhs - mat * u
        w.data = pre * d
        s.data = w
        wdn = w.InnerProduct(d, conjugate=conjugate)
        if wdn == 0:
            wdn = 1.
        err0 = sqrt(abs(wdn))

        self.errors = [err0]
        if wdn==err0:
            return u
        lwstart = log(err0)
        errstop = err0 * tol
        if self.abstol is not None:
            errstop = max(errstop, self.abstol)
        logerrstop = log(errstop)

        for it in range(maxsteps):
            self.iterations = it+1
            w.data = mat * s
            wd = wdn
            as_s = s.InnerProduct(w, conjugate=conjugate)        
            if as_s == 0: break
            alpha = wd / as_s
            u.data += alpha * s
            d.data += (-alpha) * w

            w.data = pre*d

            wdn = w.InnerProduct(d, conjugate=conjugate)
            beta = wdn / wd

            s *= beta
            s.data += w

            err = sqrt(abs(wd))
            self.errors.append(err)
            if self.printing:
                print("iteration " + str(it) + " error = " + str(err))
            if callback is not None:
                callback(it,err)
            _SetThreadPercentage(100.*max(it/maxsteps, (log(err)-lwstart)/(logerrstop - lwstart)))
            if err < errstop: break
        else:
            if self.printing:
                print("CG did not converge to tol")
        if old_status[0] != "idle":
            _PushStatus(old_status[0])
            _SetThreadPercentage(old_status[1])
            
        


def CG(mat, rhs, pre=None, sol=None, tol=1e-12, maxsteps = 100, printrates = True, initialize = True, conjugate=False, callback=None):
    """preconditioned conjugate gradient method


    Parameters
    ----------

    mat : Matrix
      The left hand side of the equation to solve. The matrix has to be spd or hermitsch.

    rhs : Vector
      The right hand side of the equation.

    pre : Preconditioner
      If provided the preconditioner is used.

    sol : Vector
      Start vector for CG method, if initialize is set False. Gets overwritten by the solution vector. If sol = None then a new vector is created.

    tol : double
      Tolerance of the residuum. CG stops if tolerance is reached.

    maxsteps : int
      Number of maximal steps for CG. If the maximal number is reached before the tolerance is reached CG stops.

    printrates : bool
      If set to True then the error of the iterations is displayed.

    initialize : bool
      If set to True then the initial guess for the CG method is set to zero. Otherwise the values of the vector sol, if provided, is used.

    conjugate : bool
      If set to True, then the complex inner product is used.


    Returns
    -------
    (vector)
      Solution vector of the CG method.

    """
    solver = CGSolver(mat=mat, pre=pre, conjugate=conjugate, tol=tol, maxsteps=maxsteps,
                      callback=callback, printing=printrates)
    solver.Solve(rhs=rhs, sol=sol, initialize=initialize)
    return solver.sol



@TimeFunction
def QMR(mat, rhs, fdofs, pre1=None, pre2=None, sol=None, maxsteps = 100, printrates = True, initialize = True, ep = 1.0, tol = 1e-7):
    """Quasi Minimal Residuum method


    Parameters
    ----------

    mat : Matrix
      The left hand side of the equation to solve

    rhs : Vector
      The right hand side of the equation.

    fdofs : BitArray
      BitArray of free degrees of freedoms.

    pre1 : Preconditioner
      First preconditioner if provided

    pre2 : Preconditioner
      Second preconditioner if provided

    sol : Vector
      Start vector for QMR method, if initialize is set False. Gets overwritten by the solution vector. If sol = None then a new vector is created.

    maxsteps : int
      Number of maximal steps for QMR. If the maximal number is reached before the tolerance is reached QMR stops.

    printrates : bool
      If set to True then the error of the iterations is displayed.

    initialize : bool
      If set to True then the initial guess for the QMR method is set to zero. Otherwise the values of the vector sol, if provided, is used.

    ep : double
      Start epsilon.

    tol : double
      Tolerance of the residuum. QMR stops if tolerance is reached.


    Returns
    -------
    (vector)
      Solution vector of the QMR method.

    """
    u = sol if sol else rhs.CreateVector()
    
    r = rhs.CreateVector()
    v = rhs.CreateVector()  
    v_tld = rhs.CreateVector()
    w = rhs.CreateVector()
    w_tld = rhs.CreateVector()
    y = rhs.CreateVector()
    y_tld = rhs.CreateVector()
    z = rhs.CreateVector()
    z_tld = rhs.CreateVector()  
    p = rhs.CreateVector()
    p_tld = rhs.CreateVector()
    q = rhs.CreateVector()  
    d = rhs.CreateVector()
    s = rhs.CreateVector()

    if (initialize): u[:] = 0.0

    r.data = rhs - mat * u
    v_tld.data = r
    y.data = pre1 * v_tld if pre1 else v_tld
    
    rho = InnerProduct(y,y)
    rho = sqrt(rho)
    
    w_tld.data = r 
    z.data = pre2.T * w_tld if pre2 else w_tld
    
    xi = InnerProduct(z,z)
    xi = sqrt(xi)
    
    gamma = 1.0
    eta = -1.0
    theta = 0.0

    
    for i in range(1,maxsteps+1):
        if (rho == 0.0):
            print('Breakdown in rho')
            return
        if (xi == 0.0):
            print('Breakdown in xi')
            return
        v.data = (1.0/rho) * v_tld
        y.data = (1.0/rho) * y

        w.data = (1.0/xi) * w_tld
        z.data = (1.0/xi) * z

        delta = InnerProduct(z,y)
        if (delta == 0.0):
            print('Breakdown in delta')
            return
        
        y_tld.data = pre2 * y if pre2 else y
        z_tld.data = pre1.T * z if pre1 else z

        if (i > 1):
            p.data = (-xi*delta / ep) * p
            p.data += y_tld

            q.data = (-rho * delta / ep) * q
            q.data += z_tld
        else:
            p.data = y_tld
            q.data = z_tld
        
        p_tld.data = mat * p
        ep = InnerProduct(q, p_tld)
        if (ep == 0.0):
            print('Breakdown in epsilon')
            return

        beta = ep/delta
        if (beta == 0.0):
            print('Breakdown in beta')
            return
        
        v_tld.data = p_tld - beta * v;      

        y.data = pre1 * v_tld if pre1 else v_tld
        
        rho_1 = rho
        rho = InnerProduct(y,y)
        rho = sqrt(rho)     
        
        
        w_tld.data = mat.T * q
        w_tld.data -= beta * w
        
        z.data = pre2.T * w_tld if pre2 else w_tld      
        
        xi = InnerProduct(z,z)
        xi = sqrt(xi)

        gamma_1 = gamma
        theta_1 = theta

        theta = rho/(gamma_1 * abs(beta)) 
        gamma = 1.0 / sqrt(1.0 + theta * theta)
        if (gamma == 0.0):
            print('Breakdown in gamma')
            return
    
        eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);        

        if (i > 1):
            d.data = (theta_1 * theta_1 * gamma * gamma) * d
            d.data += eta * p

            s.data = (theta_1 * theta_1 * gamma * gamma) * s    
            s.data += eta * p_tld

        else: 
            d.data = eta * p
            s.data = eta * p_tld
        
        
        u.data += d
        r.data -= s

        #Projected residuum: Better terminating condition necessary?
        ResNorm = sqrt( np.dot(r.FV().NumPy()[fdofs],r.FV().NumPy()[fdofs]))
        #ResNorm = sqrt(InnerProduct(r,r))

        if (printrates):
            print ("it = ", i, " err = ", ResNorm)
            
        if (ResNorm <= tol):
            break
    else:
        print("Warning: QMR did not converge to TOL")

    return u




#Source: Michael Kolmbauer https://www.numa.uni-linz.ac.at/Teaching/PhD/Finished/kolmbauer-diss.pdf
@TimeFunction
def MinRes(mat, rhs, pre=None, sol=None, maxsteps = 100, printrates = True, initialize = True, tol = 1e-7):
    """Minimal Residuum method


    Parameters
    ----------

    mat : Matrix
      The left hand side of the equation to solve

    rhs : Vector
      The right hand side of the equation.

    pre : Preconditioner
      If provided the preconditioner is used.

    sol : Vector
      Start vector for MinRes method, if initialize is set False. Gets overwritten by the solution vector. If sol = None then a new vector is created.

    maxsteps : int
      Number of maximal steps for MinRes. If the maximal number is reached before the tolerance is reached MinRes stops.

    printrates : bool
      If set to True then the error of the iterations is displayed.

    initialize : bool
      If set to True then the initial guess for the MinRes method is set to zero. Otherwise the values of the vector sol, if prevented, is used.

    tol : double
      Tolerance of the residuum. MinRes stops if tolerance is reached.


    Returns
    -------
    (vector)
      Solution vector of the MinRes method.

    """
    u = sol if sol else rhs.CreateVector()

    v_new = rhs.CreateVector()
    v = rhs.CreateVector()  
    v_old = rhs.CreateVector()
    w_new = rhs.CreateVector()
    w = rhs.CreateVector()
    w_old = rhs.CreateVector()
    z_new = rhs.CreateVector()
    z = rhs.CreateVector()
    mz = rhs.CreateVector()


    if (initialize):
        u[:] = 0.0
        v.data = rhs
    else:
        v.data = rhs - mat * u
        
    z.data = pre * v if pre else v
    
    #First Step
    gamma = sqrt(InnerProduct(z,v))
    gamma_new = 0
    z.data = 1/gamma * z 
    v.data = 1/gamma * v    
    
    ResNorm = gamma      
    ResNorm_old = gamma  
    
    if (printrates):
        print("it = ", 0, " err = ", ResNorm)       
        
    eta_old = gamma
    c_old = 1
    c = 1
    s_new = 0
    s = 0
    s_old = 0

    v_old[:] = 0.0
    w_old[:] = 0.0
    w[:] = 0.0
    
    k = 1
    while (k < maxsteps+1 and ResNorm > tol):
        mz.data = mat*z
        delta = InnerProduct(mz,z)
        v_new.data = mz - delta*v - gamma * v_old
        
        z_new.data = pre * v_new if pre else v_new

        gamma_new = sqrt(InnerProduct(z_new, v_new))
        z_new *= 1/gamma_new
        v_new *= 1/gamma_new
        
        alpha0 = c*delta - c_old*s*gamma    
        alpha1 = sqrt(alpha0*alpha0 + gamma_new*gamma_new) #**
        alpha2 = s*delta + c_old*c*gamma
        alpha3 = s_old * gamma

        c_new = alpha0/alpha1
        s_new = gamma_new/alpha1

        w_new.data = z - alpha3*w_old - alpha2*w
        w_new.data = 1/alpha1 * w_new   

        u.data += c_new*eta_old * w_new
        eta = -s_new * eta_old

        
        #update of residuum
        ResNorm = abs(s_new) * ResNorm_old
        if (printrates):        
            print("it = ", k, " err = ", ResNorm)   
        if ResNorm < tol: break
        k += 1

        # shift vectors by renaming
        v_old, v, v_new = v, v_new, v_old
        w_old, w, w_new = w, w_new, w_old
        z, z_new = z_new, z
        
        eta_old = eta
        
        s_old = s
        s = s_new

        c_old = c
        c = c_new

        gamma = gamma_new
        ResNorm_old = ResNorm
    else:
        print("Warning: MinRes did not converge to TOL")

    return u


@TimeFunction
def PreconditionedRichardson(a, rhs, pre=None, freedofs=None, maxit=100, tol=1e-8, dampfactor=1.0, printing=True):
    """ Preconditioned Richardson Iteration

    Parameters
    ----------
    a : BilinearForm
      The left hand side of the equation to solve

    rhs : Vector
      The right hand side of the equation.

    pre : Preconditioner
      If provided the preconditioner is used.
    
    freedofs : BitArray
      The FreeDofs on which the Richardson iteration acts. If argument is 'None' then the FreeDofs of the underlying FESpace is used.

    maxit : int
      Number of maximal iteration for Richardson iteration. If the maximal number is reached before the tolerance is reached a warning is displayed.

    tol : double
      Tolerance of the residuum. Richardson iteration stops if residuum < tolerance*initial_residuum is reached.

    dampfactor : float
      Set the damping factor for the Richardson iteration. If it is 1 then no damping is done. Values greater than 1 are allowed.

    printing : bool
      Set if Richardson iteration should print informations about the actual iteration like the residuum. 

    Returns
    -------
    (vector)
      Solution vector of the Preconditioned Richardson iteration.

    """

    u = rhs.CreateVector()
    r = rhs.CreateVector()
    u[:] = 0

    projector = Projector(freedofs if freedofs else a.space.FreeDofs(coupling=a.condense), False)
    
    r.data = rhs # r.data = rhs - a.mat*u
    r.data -= projector*r

    it = 0
    initial_res_norm = Norm(r)
    if printing:
        print("it =", it, " ||res||_2 =", initial_res_norm) 
    
    for it in range(1, maxit+1):
        u.data += dampfactor*(pre*r if pre else r)
        r.data = rhs - a.mat*u
        r.data -= projector*r

        res_norm = Norm(r)
        if printing:
            print("it =", it, " ||res||_2 =", res_norm)  
        if res_norm < tol*initial_res_norm:
            break
    else:
        print("Warning: Preconditioned Richardson did not converge to TOL")    

    return u   


def GMRes(A, b, pre=None, freedofs=None, x=None, maxsteps = 100, tol = None, innerproduct=None,
          callback=None, restart=None, startiteration=0, printrates=True, reltol=None):
    """ Restarting preconditioned gmres solver for A*x=b. Minimizes the preconditioned
residuum pre*(b-A*x)

Parameters
----------

A : BaseMatrix
  The left hand side of the linear system.

b : BaseVector
  The right hand side of the linear system.

pre : BaseMatrix = None
  The preconditioner for the system. If no preconditioner is given, the freedofs
  of the system must be given.

freedofs : BitArray = None
  Freedofs to solve on, only necessary if no preconditioner is given.

x : BaseVector = None
  Startvector, if given it will be modified in the routine and returned. Will be created
  if not given.

maxsteps : int = 100
  Maximum iteration steps.

tol : float = 1e-7
  Tolerance to be computed to. Gmres breaks if norm(pre*(b-A*x)) < tol.

innerproduct : function = None
  Innerproduct to be used in iteration, all orthogonalizations/norms are computed with
  respect to that inner product.

callback : function = None
  If given, this function is called with the solution vector x in each step. Only for debugging
  purposes, since it requires the overhead of computing x in every step.

restart : int = None
  If given, gmres is restarted with the current solution x every 'restart' steps.

startiteration : int = 0
  Internal value to count total number of iterations in restarted setup, no user input required
  here.

printrates : bool = True
  Print norm of preconditioned residual in each step.
"""

    if not innerproduct:
        innerproduct = lambda x,y: y.InnerProduct(x, conjugate=True)
        norm = Norm
    else:
        norm = lambda x: sqrt(innerproduct(x,x).real)
    is_complex = b.is_complex
    if not pre:
        assert freedofs is not None
        pre = Projector(freedofs, True)
    n = len(b)
    m = maxsteps
    if not x:
        x = b.CreateVector()
        x[:] = 0

    if callback:
        xstart = x.CreateVector()
        xstart.data = x
    else:
        xstart = None
    sn = Vector(m, is_complex)
    cs = Vector(m, is_complex)
    sn[:] = 0
    cs[:] = 0

    r = b.CreateVector()
    tmp = b.CreateVector()
    tmp.data = b - A * x
    r.data = pre * tmp

    Q = []
    H = []
    Q.append(b.CreateVector())
    r_norm = norm(r)
    if printrates:
        print("Step 0, error = ", r_norm)
    if reltol is not None and r_norm != 0:
        rtol = reltol * abs(r_norm)
        tol = rtol if tol is None else max(rtol, tol)
    if tol is None:
        tol = 1e-7
    if abs(r_norm) < tol:
        return x
    Q[0].data = 1./r_norm * r
    beta = Vector(m+1, is_complex)
    beta[:] = 0
    beta[0] = r_norm

    def arnoldi(A,Q,k):
        q = b.CreateVector()
        tmp.data = A * Q[k]
        q.data = pre * tmp
        h = Vector(m+1, is_complex)
        h[:] = 0
        for i in range(k+1):
            h[i] = innerproduct(Q[i],q)
            q.data += (-1)* h[i] * Q[i]
        h[k+1] = norm(q)
        if abs(h[k+1]) < 1e-12:
            return h, None
        q *= 1./h[k+1].real
        return h, q

    def givens_rotation(v1,v2):
        if v2 == 0:
            return 1,0
        elif v1 == 0:
            return 0,v2/abs(v2)
        else:
            t = sqrt((v1.conjugate()*v1+v2.conjugate()*v2).real)
            cs = abs(v1)/t
            sn = v1/abs(v1) * v2.conjugate()/t
            return cs,sn

    def apply_givens_rotation(h, cs, sn, k):
        for i in range(k):
            temp = cs[i] * h[i] + sn[i] * h[i+1]
            h[i+1] = -sn[i].conjugate() * h[i] + cs[i].conjugate() * h[i+1]
            h[i] = temp
        cs[k], sn[k] = givens_rotation(h[k], h[k+1])
        h[k] = cs[k] * h[k] + sn[k] * h[k+1]
        h[k+1] = 0

    def calcSolution(k):
        mat = Matrix(k+1,k+1, is_complex)
        for i in range(k+1):
            mat[:,i] = H[i][:k+1]
        rs = Vector(k+1, is_complex)
        rs[:] = beta[:k+1]
        y = mat.I * rs
        if xstart:
            x.data = xstart
        for i in range(k+1):
            x.data += y[i] * Q[i]

    for k in range(m):
        startiteration += 1
        h,q = arnoldi(A,Q,k)
        H.append(h)
        if q is None:
            break
        Q.append(q)
        apply_givens_rotation(h, cs, sn, k)
        beta[k+1] = -sn[k].conjugate() * beta[k]
        beta[k] = cs[k] * beta[k]
        error = abs(beta[k+1])
        if printrates:
            print("Step", startiteration, ", error = ", error)
        if callback:
            calcSolution(k)
            callback(x)
        if error < tol:
            break
        if restart and k+1 == restart and not (restart == maxsteps):
            calcSolution(k)
            del Q
            return GMRes(A, b, freedofs=freedofs, pre=pre, x=x, maxsteps=maxsteps-restart, callback=callback,
                         tol=tol, innerproduct=innerproduct,
                         restart=restart, startiteration=startiteration, printrates=printrates)
    calcSolution(k)
    return x
