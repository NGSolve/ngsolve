def BVP(bf, lf, gf, pre=None, maxsteps=200, tol=1e-8, print=True, inverse="umfpack"):
    """
    Solve a linear boundary value problem
    """
    from ngsolve import Projector
    from ngsolve.krylovspace import CG
    
    r = lf.vec.CreateVector()
    r.data = lf.vec
    
    if bf.condense:
        r.data += bf.harmonic_extension_trans * r

        # zero local dofs
        innerbits = gf.space.FreeDofs(False) & ~gf.space.FreeDofs(True)
        Projector(innerbits, False).Project(gf.vec)
        
    if pre:
        CG(mat = bf.mat, rhs = r, pre=pre, sol=gf.vec, tol=tol, initialize=False, printrates=print)
    else:
        inv = bf.mat.Inverse(gf.space.FreeDofs(bf.condense), inverse=inverse)
        r.data -= bf.mat * gf.vec
        gf.vec.data += inv * r

    if bf.condense:
        gf.vec.data += bf.harmonic_extension * gf.vec
        gf.vec.data += bf.inner_solve * r
            
