
from ngsolve import Projector, Norm, TimeFunction, BaseMatrix, Preconditioner, InnerProduct, \
    Norm, sqrt, Vector, Matrix, BaseVector, BlockVector, BitArray
from typing import Optional, Callable, Union
import logging
from netgen.libngpy._meshing import _PushStatus, _GetStatus, _SetThreadPercentage
from math import log
import os


if os.name == "nt":
    _clear_line_command = ""
else:
    _clear_line_command = "\33[2K"

linear_solver_param_doc = """
mat : BaseMatrix
  The left hand side of the equation to solve.

pre : Preconditioner, BaseMatrix = None
  If provided, the preconditioner for the system.

freedofs : BitArray = None
  If no preconditioner is provided, the BitArray of the FESpace freedofs must be given.

tol : double = 1e-12
  Relative tolerance for the residuum reduction.

maxiter : int = 100
  Maximum number of iterations, if reached solver will emit a warning.

callback : Callable[[int, float], None] = None
  Callback function that is called with iteration number and residual in each iteration step.

callback_sol : Callable[[BaseVector], None] = None
  Callback function that is called with solution x_k in each iteration step.

printrates : bool = False
  Print iterations to stdout. One can give a string to be passed as an `end`
  argument to the print function, for example:
  >>> printrates="\r"
  will call
  >>> print("iteration = 1, residual = 1e-3", end="\r")
  if "\r" is passed, a final output will also be printed.

plotrates : bool = False
  matplotlib plot of errors (residuals)
"""

class LinearSolver(BaseMatrix):
    """Base class for linear solvers.
""" + linear_solver_param_doc
    name = "LinearSolver"
    def __init__(self, mat : BaseMatrix,
                 pre : Optional[Preconditioner] = None,
                 freedofs : Optional[BitArray] = None,
                 tol : float = None,
                 maxiter : int = 100,
                 atol : float = None,
                 callback : Optional[Callable[[int, float], None]] = None,
                 callback_sol : Optional[Callable[[BaseVector], None]] = None,
                 printrates : bool = False,
                 plotrates : bool = False):
        super().__init__()
        if atol is None and tol is None:
            tol = 1e-12
        self.mat = mat
        assert (freedofs is None) != (pre is None) # either pre or freedofs must be given
        self.pre = pre if pre else Projector(freedofs, True)
        self.tol = tol
        self.atol = atol
        self.maxiter = maxiter
        self.callback = callback
        self.callback_sol = callback_sol
        self.printrates = printrates
        self.plotrates = plotrates
        self.residuals = []
        self.iterations = 0

    @TimeFunction
    def Solve(self, rhs : BaseVector, sol : Optional[BaseVector] = None,
              initialize : bool = True) -> BaseVector:
        self.iterations = 0
        self.residuals = []
        old_status = _GetStatus()
        _PushStatus(self.name + " Solve")
        _SetThreadPercentage(0)
        if sol is None:
            sol = rhs.CreateVector()
            initialize = True
        if initialize:
            sol[:] = 0
        self.sol = sol
        self._SolveImpl(rhs=rhs, sol=sol)
        if old_status[0] != "idle":
            _PushStatus(old_status[0])
            _SetThreadPercentage(old_status[1])
        return sol

    def Height(self) -> int:
        return self.mat.width

    def Width(self) -> int:
        return self.mat.height

    def CreateVector(self,col):
        return self.mat.CreateVector(not col)
    
    def IsComplex(self) -> bool:
        return self.mat.IsComplex()

    def Mult(self, x : BaseVector, y : BaseVector) -> None:
        self.Solve(rhs=x, sol=y, initialize=True)

    def Update(self):
        if hasattr(self.pre, "Update"):
            self.pre.Update()

    def CheckResidual(self, residual):
        self.iterations += 1
        self.residuals.append(residual)
        if len(self.residuals) == 1:
            if self.tol is None:
                self._final_residual = self.atol
            else:
                self._final_residual = residual * self.tol
                if self.atol is not None:
                    self._final_residual = max(self._final_residual, self.atol)
        else:
            if self.callback is not None:
                self.callback(self.iterations, residual)
            if self.callback_sol is not None:
                self.callback_sol(self.sol)
        if self.residuals[0] != 0:
            logerrstop = log(self._final_residual)
            logerrfirst = log(self.residuals[0])
            if residual == 0:
                _SetThreadPercentage(100)
            else:
                _SetThreadPercentage(100.*max(self.iterations/self.maxiter,
                                              (log(residual)-logerrfirst)/(logerrstop - logerrfirst)))
        if self.printrates:
            print("{}{} iteration {}, residual = {}     ".format(_clear_line_command, self.name, self.iterations, residual), end="\n" if isinstance(self.printrates, bool) else self.printrates)
            if self.iterations == self.maxiter and residual > self._final_residual:
                print("{}WARNING: {} did not converge to TOL".format(_clear_line_command, self.name))
        is_converged = self.iterations >= self.maxiter or residual <= self._final_residual
        if is_converged and self.printrates == "\r":
            print("{}{} {}converged in {} iterations to residual {}".format(_clear_line_command, self.name, "NOT " if residual >= self._final_residual else "", self.iterations, residual))

        if self.plotrates:
            if self.iterations==1:
                import matplotlib.pyplot as plt
                from IPython.display import display, clear_output
                fig, ax = plt.subplots()
                self.plt = plt
                self.ax = ax
                self.fig = fig
                self.its = []
                self.ress = []
                self.clear_output=clear_output
                self.display=display
                plt.ioff()
                plt.show()
            self.its.append(self.iterations)
            self.ress.append(residual)
            # update_plot(plt, ax, self.its, self.ress)
            self.ax.clear()
            self.ax.semilogy(self.its, self.ress, label='error')
            self.ax.set_xlabel('iteration')
            self.ax.set_ylabel('error')
            self.ax.set_title('CG Solver Convergence')
            self.ax.legend()
            self.plt.draw()
            self.clear_output(wait=True)
            self.display(self.fig)
            
            
        return is_converged

class CGSolver(LinearSolver):
    """Preconditioned conjugate gradient method

    Parameters
    ----------

""" + linear_solver_param_doc + """

conjugate : bool = False
  If set to True, then the complex inner product is used, else a pseudo inner product that makes CG work with complex symmetric matrices.
"""
    name = "CG"

    def __init__(self, *args,
                 conjugate : bool = False,
                 abstol : float = None,
                 maxsteps : int = None,
                 printing : bool = False,
                 **kwargs):
        if printing:
            print("WARNING: printing is deprecated, use printrates instead!")
            kwargs["printrates"] = printing
        if abstol is not None:
            print("WARNING: abstol is deprecated, use atol instead!")
            kwargs["abstol"] = abstol
        if maxsteps is not None:
            print("WARNING: maxsteps is deprecated, use maxiter instead!")
            kwargs["maxiter"] = maxsteps
        super().__init__(*args, **kwargs)
        self.conjugate = conjugate

    # for backward compatibility
    @property
    def errors(self):
        return self.residuals

    def _SolveImpl(self, rhs : BaseVector, sol : BaseVector):
        d, w, s = [sol.CreateVector() for i in range(3)]
        conjugate = self.conjugate
        d.data = rhs - self.mat * sol
        w.data = self.pre * d
        s.data = w
        wdn = w.InnerProduct(d, conjugate=conjugate)
        if self.CheckResidual(sqrt(abs(wdn))):
            return

        while True:
            w.data = self.mat * s
            wd = wdn
            as_s = s.InnerProduct(w, conjugate=conjugate)        
            if as_s == 0 or wd == 0: break
            alpha = wd / as_s
            sol.data += alpha * s
            d.data += (-alpha) * w

            w.data = self.pre * d

            wdn = w.InnerProduct(d, conjugate=conjugate)
            if self.CheckResidual(sqrt(abs(wdn))):
                return

            beta = wdn / wd
            s *= beta
            s.data += w

        
def CG(mat, rhs, pre=None, sol=None, tol=1e-12, maxsteps = 100, printrates = True, plotrates = False, initialize = True, conjugate=False, callback=None, **kwargs):
    """preconditioned conjugate gradient method


    Parameters
    ----------

    mat : Matrix
      The left hand side of the equation to solve. The matrix has to be spd o hermitsch.

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

    plotrates : bool
      If set to True then the error of the iterations is plotted.

    initialize : bool
      If set to True then the initial guess for the CG method is set to zero. Otherwise the values of the vector sol, if provided, is used.

    conjugate : bool
      If set to True, then the complex inner product is used.


    Returns
    -------
    (vector)
      Solution vector of the CG method.

    """
    solver = CGSolver(mat=mat, pre=pre, conjugate=conjugate, tol=tol, maxiter=maxsteps,
                      callback=callback, printrates=printrates, plotrates=plotrates, **kwargs)
    solver.Solve(rhs=rhs, sol=sol, initialize=initialize)
    return solver.sol

class QMRSolver(LinearSolver):
    """Quasi Minimal Residuum method

    Parameters
    ----------

""" + linear_solver_param_doc + """

pre2 : Preconditioner = None
  Second preconditioner, if provided.

ep : double
  Start epsilon.
"""

    name = "QMR"

    def __init__(self, *args, pre2 : Preconditioner = None,
                 ep : float = 1., **kwargs):
        super().__init__(*args, **kwargs)
        self.pre2 = pre2
        self.ep = ep

    def _SolveImpl(self, rhs : BaseVector, sol : BaseVector):
        u, mat, ep, pre1, pre2 = sol, self.mat, self.ep, self.pre, self.pre2
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

        r.data = rhs - mat * u
        v_tld.data = r
        y.data = pre1 * v_tld

        rho = InnerProduct(y,y)
        rho = sqrt(rho)

        w_tld.data = r
        z.data = pre2.T * w_tld if pre2 else w_tld

        xi = InnerProduct(z,z)
        xi = sqrt(xi)

        gamma = 1.0
        eta = -1.0
        theta = 0.0

        for i in range(1,self.maxiter+1):
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
            z_tld.data = pre1.T * z

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

            y.data = pre1 * v_tld

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
            v.data = self.pre * r
            ResNorm = sqrt(InnerProduct(r,v))
            # ResNorm = sqrt( np.dot(r.FV().NumPy()[fdofs],r.FV().NumPy()[fdofs]))
            #ResNorm = sqrt(InnerProduct(r,r))
            if self.CheckResidual(ResNorm):
                return

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
    # backwards compatibility, but freedofs are not needed then.
    if pre1 is not None:
        fdofs = None
    return QMRSolver(mat=mat, freedofs=fdofs, pre=pre1,
                     pre2=pre2, maxiter=maxsteps,
                     printrates=printrates, ep=ep,
                     tol=tol).Solve(rhs=rhs, sol=sol, initialize=initialize)




#Source: Michael Kolmbauer https://www.numa.uni-linz.ac.at/Teaching/PhD/Finished/kolmbauer-diss.pdf
class MinResSolver(LinearSolver):
    """Minimal Residuum method

    Parameters
    ----------
""" + linear_solver_param_doc
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _SolveImpl(self, rhs: BaseVector, sol : BaseVector):
        pre, mat, u = self.pre, self.mat, sol
        v_new = rhs.CreateVector()
        v = rhs.CreateVector()
        v_old = rhs.CreateVector()
        w_new = rhs.CreateVector()
        w = rhs.CreateVector()
        w_old = rhs.CreateVector()
        z_new = rhs.CreateVector()
        z = rhs.CreateVector()
        mz = rhs.CreateVector()

        v.data = rhs - mat * u

        z.data = pre * v

        #First Step
        gamma = sqrt(InnerProduct(z,v))
        gamma_new = 0
        z.data = 1/gamma * z
        v.data = 1/gamma * v

        ResNorm = gamma
        ResNorm_old = gamma

        if self.CheckResidual(ResNorm):
            return

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
        while True:
            mz.data = mat*z
            delta = InnerProduct(mz,z)
            v_new.data = mz - delta*v - gamma * v_old

            z_new.data = pre * v_new

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
            if self.CheckResidual(ResNorm):
                return
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
    return MinResSolver(mat=mat, pre=pre, maxiter=maxsteps,
                        printrates=printrates,
                        tol=tol).Solve(rhs=rhs, sol=sol,
                                       initialize=initialize)


class RichardsonSolver(LinearSolver):
    """ Preconditioned Richardson Iteration

Parameters
----------
""" + linear_solver_param_doc + """

dampfactor : float = 1.
  Set the damping factor for the Richardson iteration. If it is 1 then no damping is done. Values greater than 1 are allowed.
"""
    name = "Richardson"
    def __init__(self, *args, dampfactor = 1., **kwargs):
        super().__init__(*args, **kwargs)
        self.dampfactor = dampfactor

    def _SolveImpl(self, rhs : BaseVector, sol : BaseVector):
        r = rhs.CreateVector()
        d = sol.CreateVector()
        r.data = rhs - self.mat*sol
        d.data = self.pre * r
        res_norm = abs(InnerProduct(d,r))
        if self.CheckResidual(res_norm):
            return

        while True:
            sol.data += self.dampfactor * d
            r.data = rhs - self.mat * sol
            d.data = self.pre * r

            res_norm = abs(InnerProduct(d,r))
            if self.CheckResidual(res_norm):
                return


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

class GMResSolver(LinearSolver):
    """Preconditioned GMRes solver. Minimizes the preconditioned residuum pre * (b-A*x)

Parameters
----------

""" + linear_solver_param_doc + """

innerproduct : Callable[[BaseVector, BaseVector], Union[float, complex]] = None
  Innerproduct to be used in iteration, all orthogonalizations/norms are computed with respect to that inner product.

restart : int = None
  If given, GMRes is restarted with the current solution x every 'restart' steps.
"""
    name = "GMRes"

    def __init__(self, *args,
                 innerproduct : Optional[Callable[[BaseVector, BaseVector],
                                                  Union[float, complex]]] = None,
                 restart : Optional[int] = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        if innerproduct is not None:
            self.innerproduct = innerproduct
            self.norm = lambda x: sqrt(innerproduct(x,x).real)
            self.restart = restart
        else:
            self.innerproduct = lambda x, y: y.InnerProduct(x, conjugate=True)
            self.norm = Norm
            self.restart = restart

    def _SolveImpl(self, rhs : BaseVector, sol : BaseVector):
        is_complex = rhs.is_complex
        A, pre, innerproduct, norm = self.mat, self.pre, self.innerproduct, self.norm
        n = len(rhs)
        m = self.maxiter
        sn = Vector(m, is_complex)
        cs = Vector(m, is_complex)
        sn[:] = 0
        cs[:] = 0
        if self.callback_sol is not None:
            sol_start = sol.CreateVector()
            sol_start.data = sol
        r = rhs.CreateVector()
        tmp = rhs.CreateVector()
        tmp.data = rhs - A * sol
        r.data = pre * tmp
        Q = []
        H = []
        Q.append(rhs.CreateVector())
        r_norm = norm(r)
        if self.CheckResidual(abs(r_norm)):
            return sol
        Q[0].data = 1./r_norm * r
        beta = Vector(m+1, is_complex)
        beta[:] = 0
        beta[0] = r_norm

        def arnoldi(A,Q,k):
            q = rhs.CreateVector()
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
            # if callback_sol is set we need to recompute solution in every step
            if self.callback_sol is not None:
                sol.data = sol_start
            mat = Matrix(k+1,k+1, is_complex)
            for i in range(k+1):
                mat[:,i] = H[i][:k+1]
            rs = Vector(k+1, is_complex)
            rs[:] = beta[:k+1]
            y = mat.I * rs
            for i in range(k+1):
                sol.data += y[i] * Q[i]

        for k in range(m):
            h,q = arnoldi(A,Q,k)
            H.append(h)
            if q is None:
                break
            Q.append(q)
            apply_givens_rotation(h, cs, sn, k)
            beta[k+1] = -sn[k].conjugate() * beta[k]
            beta[k] = cs[k] * beta[k]
            error = abs(beta[k+1])
            if self.callback_sol is not None:
                calcSolution(k)
            if self.CheckResidual(error):
                break
            if self.restart is not None and (k+1 == self.restart and not (self.restart == self.maxiter)):
                calcSolution(k)
                del Q
                restarted_solver = GMResSolver(mat=self.mat,
                                               pre=self.pre,
                                               tol=0,
                                               atol=self._final_residual,
                                               callback=self.callback,
                                               callback_sol=self.callback_sol,
                                               maxiter=self.maxiter,
                                               restart=self.restart,
                                               printrates=self.printrates)
                restarted_solver.iterations = self.iterations
                sol = restarted_solver.Solve(rhs = rhs, sol = sol, initialize=False)
                self.residuals += restarted_solver.residuals
                self.iterations = restarted_solver.iterations
                return sol
        calcSolution(k)
        return sol

def GMRes(A, b, pre=None, freedofs=None, x=None, maxsteps = 100, tol = None, innerproduct=None,
          callback=None, restart=None, startiteration=0, printrates=True, reltol=None):
    """Restarting preconditioned gmres solver for A*x=b. Minimizes the preconditioned residuum pre*(b-A*x).

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

innerproduct : function = None
  Innerproduct to be used in iteration, all orthogonalizations/norms are computed with
  respect to that inner product.

callback : function = None
  If given, this function is called with the solution vector x in each step. Only for debugging

restart : int = None
  If given, gmres is restarted with the current solution x every 'restart' steps.

startiteration : int = 0
  Internal value to count total number of iterations in restarted setup, no user input required
  here.

printrates : bool = True
  Print norm of preconditioned residual in each step.
"""
    solver = GMResSolver(mat=A, pre=pre, freedofs=freedofs,
                         maxiter=maxsteps, tol=reltol, atol=tol,
                         innerproduct=innerproduct,
                         callback_sol=callback, restart=restart,
                         printrates=printrates)
    return solver.Solve(rhs=b, sol=x)




from ngsolve.la import EigenValues_Preconditioner
def BramblePasciakCG(A, B, C, f, g, preA, preS, maxit=1000, tol=1e-8, \
                         printrates=False):

    printeol = "\n"
    if isinstance(printrates, str):
        printeol = printrates
        printrates = True

    lam = EigenValues_Preconditioner(A,preA)
    if printrates==True:
        print ("lammin/lammax = ", lam[0], '/', lam[-1])
    preA = 1.2/lam[0]*preA   # scaling

    
    x = BlockVector([f.CreateVector(), g.CreateVector()])
    w,r,p,ap = [x.CreateVector() for i in range(4)]
    
    # r.data = b
    # p.data = pre*r
    pru = (preA * f).Evaluate()
    r[0].data = A*pru - f
    r[1].data = B*pru - g
    p[0].data = pru    
    p[1].data = preS*r[1]
    
    wrn = InnerProduct(r,p)
    err0 = sqrt(wrn)
        
    x[:] = 0
    for it in range(maxit):
        # ap.data = A * p
        hv = (A * p[0] + B.T * p[1]).Evaluate()
        papu = (preA * hv).Evaluate()
        ap[0].data = A * papu - hv
        ap[1].data = B * (papu - p[0])
        if C is not None:
            ap[1].data += C * p[1]
        
        pap = InnerProduct(p, ap)
        wr = wrn
        alpha = wr / pap
        
        x += alpha * p
        r -= alpha * ap
        pru -= alpha * papu

        # w.data = pre*r
        w[0].data = pru
        w[1].data = preS * r[1]
        
        wrn = InnerProduct(w, r)
        err = sqrt(wrn)
        if printrates==True:
            print ("Iteration",it,"err=",err,"    ",end=printeol)
        if err < tol * err0: break
            
        beta = wrn / wr
        
        p *= beta
        p.data += w 
    
    return x[0], x[1]



def update_plot(plt, ax, its, ress):
    # its.append(it)
    # ress.append(res)
    ax.clear()
    ax.semilogy(its, ress, label='error')
    ax.set_xlabel('iteration')
    ax.set_ylabel('error')
    ax.set_title('CG Solver Convergence')
    ax.legend()
    plt.draw()
    clear_output(wait=True)
    display(fig)
