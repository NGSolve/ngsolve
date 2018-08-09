
def Newton(a, u, freedofs=None, maxit=100, maxerr=1e-11, inverse="umfpack", el_int=False, dampfactor=1, printing=True):
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

    el_int : bool
      Flag if eliminate internal flag is used. If this is set to True then Newton uses static condensation to invert the matrix. If freedofs is not None, the user has to take care that the local dofs are set to zero in the freedofs array.

    dampfactor : float
      Set the damping factor for Newton's method. If dampfactor is 1 then no damping is done. If value is < 1 then the damping is done by the formula 'min(1,dampfactor*numit)' for the correction, where 'numit' denotes the Newton iteration.

    printing : bool
      Set if Newton's method should print informations about the actual iteration like the error. 

    Returns
    -------
    (int, int)
      List of two integers. The first one is 0 if Newton's method did converge, -1 otherwise. The second one gives the number of Newton iterations needed.

    """
    w = u.vec.CreateVector()
    r = u.vec.CreateVector()

    err = 1
    numit = 0
    inv = None
    
    for it in range(maxit):
        numit += 1
        if printing:
            print("Newton iteration ", it)
            print ("energy = ", a.Energy(u.vec))
        a.Apply(u.vec, r)
        a.AssembleLinearization(u.vec)

        inv = a.mat.Inverse(freedofs if freedofs else u.space.FreeDofs(el_int), inverse=inverse)

        if el_int:
            r.data += a.harmonic_extension_trans * r
            w.data = inv * r
            w.data += a.harmonic_extension * w
            w.data += a.inner_solve * r
        else:
            w.data = inv * r
            
        err = sqrt(abs(InnerProduct(w,r)))
        if printing:
            print("err = ", err)

        energy = a.Energy(u.vec)
        u.vec.data -= min(1, numit*dampfactor)*w
        
        if abs(err) < maxerr: break
    else:
        print("Warning: Newton might not converge!")
        return (-1,numit)
    return (0,numit)




def NewtonMinimization(a, u, freedofs=None, maxit=100, maxerr=1e-11, inverse="umfpack", el_int=False, dampfactor=1, linesearch=False, printing=True):
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

    el_int : bool
      Flag if eliminate internal flag is used. If this is set to True then Newton uses static condensation to invert the matrix. If freedofs is not None, the user has to take care that the local dofs are set to zero in the freedofs array.

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
    w = u.vec.CreateVector()
    r = u.vec.CreateVector()
    uh = u.vec.CreateVector()

    err = 1
    numit = 0
    inv = None
    
    for it in range(maxit):
        numit += 1
        if printing:
            print("Newton iteration ", it)
            print ("energy = ", a.Energy(u.vec))
        a.Apply(u.vec, r)
        a.AssembleLinearization(u.vec)

        inv = a.mat.Inverse(freedofs if freedofs else u.space.FreeDofs(el_int), inverse=inverse)

        if el_int:
            r.data += a.harmonic_extension_trans * r
            w.data = inv * r
            w.data += a.harmonic_extension * w
            w.data += a.inner_solve * r
        else:
            w.data = inv * r
            
        err = sqrt(abs(InnerProduct(w,r)))
        if printing:
            print("err = ", err)

        energy = a.Energy(u.vec)
        uh.data = u.vec - min(1, numit*dampfactor)*w
            
        tau = min(1, numit*dampfactor)
        if linesearch:
            while a.Energy(uh) > energy+1e-15:
                tau *= 0.5
                uh.data = u.vec - tau * w
                if printing:
                    print ("tau = ", tau)
                    print ("energy uh = ", a.Energy(uh))
        u.vec.data = uh
        if abs(err) < maxerr: break
    else:
        print("Warning: Newton might not converge!")
        return (-1,numit)
    return (0,numit)



def PreconditionedRichardson(a, rhs, pre=None, freedofs=None, maxit=100, tol=1e-8, el_int=False, dampfactor=1.0, printing=True):
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

    el_int : bool
      Flag if eliminate internal flag is used. If this is set to True then Richardson iteration acts only on the non-local dofs. If freedofs is not None, the user has to take care that the local dofs are set to zero in the freedofs array.

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

    projector = Projector(freedofs if freedofs else a.space.FreeDofs(coupling=el_int), False)
    
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


