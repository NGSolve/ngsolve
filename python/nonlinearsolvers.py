from ngsolve.la import InnerProduct
from math import sqrt
from ngsolve import Projector, Norm
from .utils import TimeFunction

class NewtonSolver:
    def __init__(self, a, u, rhs=None, freedofs=None,
                 inverse="umfpack", solver=None):
        self.a, self.u, self.inverse = a, u, inverse
        self.w = u.vec.CreateVector()
        self.r = u.vec.CreateVector()
        self.rhs = rhs
        self.uh = u.vec.CreateVector()
        self.inv = None if solver is None else solver
        if solver:
            self.inverse = "given"
        else:
            self.freedofs = freedofs or u.space.FreeDofs(a.condense)

    @TimeFunction
    def Solve(self, maxit=100, maxerr=1e-11, dampfactor=1,
              printing=False, callback=None, linesearch=False,
              printenergy=False, print_wrong_direction=False):
        numit = 0
        err = 1.
        a, u, w, r,uh = self.a, self.u, self.w, self.r, self.uh
        for it in range(maxit):
            numit += 1
            if printing:
                print("Newton iteration ", it)
                if printenergy:
                    print("Energy: ", a.Energy(u.vec))

            a.AssembleLinearization(u.vec)
            a.Apply(u.vec, r)

            self._UpdateInverse()
            if self.rhs is not None:
                r.data -= self.rhs.vec
            if a.condense:
                r.data += a.harmonic_extension_trans * r
                w.data = self.inv * r
                w.data += a.harmonic_extension * w
                w.data += a.inner_solve * r
            else:
                w.data = self.inv * r

            err2 = InnerProduct(w,r)
            if print_wrong_direction:
                if err2 < 0:
                    print("wrong direction")
            err = sqrt(abs(err2))
            if printing:
                print("err = ", err)

            tau = min(1, numit*dampfactor)

            if linesearch:
                uh.data = u.vec - tau*w
                energy = a.Energy(u.vec)
                while a.Energy(uh) > energy+(max(1e-14*abs(energy),maxerr)) and tau > 1e-10:
                    tau *= 0.5
                    uh.data = u.vec - tau * w
                    if printing:
                        print ("tau = ", tau)
                        print ("energy uh = ", a.Energy(uh))
                u.vec.data = uh

            else:
                u.vec.data -= tau * w
            if callback is not None:
                callback(it, err)
            if abs(err) < maxerr: break
        else:
            print("Warning: Newton might not converge! Error = ", err)
            return (-1,numit)
        return (0,numit)

    def SetDirichlet(self, dirichletvalues):
        a, u, w, r = self.a, self.u, self.w, self.r
        a.AssembleLinearization(u.vec)
        self._UpdateInverse()
        w.data = dirichletvalues-u.vec
        r.data = a.mat * w
        w.data -= self.inv*r
        u.vec.data += w

    def _UpdateInverse(self):
        if self.inverse in ("sparsecholesky", "given") and self.inv:
            self.inv.Update()
        else:
            self.inv = self.a.mat.Inverse(self.freedofs,
                                          inverse=self.inverse)


def Newton(a, u, freedofs=None, maxit=100, maxerr=1e-11, inverse="umfpack", \
               dirichletvalues=None, dampfactor=1, printing=True, callback=None):
    """
    Newton's method for solving non-linear problems of the form A(u)=0.

    Parameters
    ----------
    a : BilinearForm
      The BilinearForm of the non-linear variational problem. It does not have to be assembled.

    u : GridFunction
      The GridFunction where the solution is saved. The values are used as initial guess for Newton's method.

    freedofs : BitArray
      The FreeDofs on which the assembled matrix is inverted. If argument is 'None' then the FreeDofs of the underlying FESpace is used.

    maxit : int
      Number of maximal iteration for Newton. If the maximal number is reached before the maximal error Newton might no converge and a warning is displayed.

    maxerr : float
      The maximal error which Newton should reach before it stops. The error is computed by the square root of the inner product of the residuum and the correction.

    inverse : string
      A string of the sparse direct solver which should be solved for inverting the assembled Newton matrix.

    dampfactor : float
      Set the damping factor for Newton's method. If dampfactor is 1 then no damping is done. If value is < 1 then the damping is done by the formula 'min(1,dampfactor*numit)' for the correction, where 'numit' denotes the Newton iteration.

    printing : bool
      Set if Newton's method should print informations about the actual iteration like the error. 

    Returns
    -------
    (int, int)
      List of two integers. The first one is 0 if Newton's method did converge, -1 otherwise. The second one gives the number of Newton iterations needed.

    """
    solver = NewtonSolver(a=a, u=u, freedofs=freedofs, inverse=inverse)
    if dirichletvalues is not None:
        solver.SetDirichlet(dirichletvalues)
    return solver.Solve(maxit=maxit, maxerr=maxerr,
                        dampfactor=dampfactor,
                        printing=printing,
                        callback=callback,
                        linesearch=False,
                        printenergy=False)


def NewtonMinimization(a, u, freedofs=None, maxit=100, maxerr=1e-11, inverse="umfpack", dampfactor=1, linesearch=False, printing=True, callback=None):
    """
    Newton's method for solving non-linear problems of the form A(u)=0 involving energy integrators.


    Parameters
    ----------
    a : BilinearForm
      The BilinearForm of the non-linear variational problem. It does not have to be assembled.

    u : GridFunction
      The GridFunction where the solution is saved. The values are used as initial guess for Newton's method.

    freedofs : BitArray
      The FreeDofs on which the assembled matrix is inverted. If argument is 'None' then the FreeDofs of the underlying FESpace is used.

    maxit : int
      Number of maximal iteration for Newton. If the maximal number is reached before the maximal error Newton might no converge and a warning is displayed.

    maxerr : float
      The maximal error which Newton should reach before it stops. The error is computed by the square root of the inner product of the residuum and the correction.

    inverse : string
      A string of the sparse direct solver which should be solved for inverting the assembled Newton matrix.

    dampfactor : float
      Set the damping factor for Newton's method. If dampfactor is 1 then no damping is done. If value is < 1 then the damping is done by the formula 'min(1,dampfactor*numit)' for the correction, where 'numit' denotes the Newton iteration.

    linesearch : bool
      If True then linesearch is used to guarantee that the energy decreases in every Newton iteration.

    printing : bool
      Set if Newton's method should print informations about the actual iteration like the error. 

    Returns
    -------
    (int, int)
      List of two integers. The first one is 0 if Newton's method did converge, -1 otherwise. The second one gives the number of Newton iterations needed.

    """
    solver = NewtonSolver(a=a, u=u, freedofs=freedofs, inverse=inverse)
    return solver.Solve(maxit=maxit, maxerr=maxerr,
                        dampfactor=dampfactor,
                        printing=printing,
                        callback=callback,
                        linesearch=linesearch,
                        printenergy=printing,
                        print_wrong_direction=True)

