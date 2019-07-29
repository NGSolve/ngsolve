
def BVP(bf, lf, gf, \
            pre=None, pre_flags={}, \
            solver=None, solver_flags={},
            maxsteps=200, tol=1e-8, print=True, inverse="umfpack",
            needsassembling=True):
    """
    Solve a linear boundary value problem A(u,v) = f(v)

Parameters
----------

bf : BilinearForm
  provides the matrix. 

lf : LinearForm
  provides the right hand side.

gf : GridFunction
  provides the solution vector

pre : Basematrix or class or string = None
  used if an iterative solver is used
  can be one of 
    * a preconditioner object
    * a preconditioner class
    * a preconditioner class name

pre_flags : dictionary = { }
  flags used to create preconditioner


    """
    from ngsolve import Projector, Preconditioner
    from ngsolve.krylovspace import CG


    if isinstance(pre,type):
        pre = pre(bf, **pre_flags)
        if not needsassembling:
            pre.Update()

    if isinstance(pre,str):
        pre = Preconditioner(bf, pre, **pre_flags)
        if not needsassembling:
            pre.Update()

        
    if needsassembling:
        bf.Assemble()
        lf.Assemble()
    
    r = lf.vec.CreateVector()
    r.data = lf.vec
    
    if bf.condense:
        r.data += bf.harmonic_extension_trans * r

        # zero local dofs
        innerbits = gf.space.FreeDofs(False) & ~gf.space.FreeDofs(True)
        Projector(innerbits, False).Project(gf.vec)

    if solver:
        inv = solver(mat=bf.mat, pre=pre, **solver_flags)
        r.data -= bf.mat * gf.vec
        gf.vec.data += inv * r
                   
    elif pre:
        CG(mat = bf.mat, rhs = r, pre=pre, sol=gf.vec, tol=tol, maxsteps=maxsteps, initialize=False, printrates=print)
    else:
        inv = bf.mat.Inverse(gf.space.FreeDofs(bf.condense), inverse=inverse)
        r.data -= bf.mat * gf.vec
        gf.vec.data += inv * r

    if bf.condense:
        gf.vec.data += bf.harmonic_extension * gf.vec
        gf.vec.data += bf.inner_solve * r
            
