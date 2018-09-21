from math import sqrt
from ngsolve import *
def SolveNonlinearMinProblem(a,gfu):
    res = gfu.vec.CreateVector()
    du  = gfu.vec.CreateVector()

    for it in range(12):
        print ("Newton iteration", it)
        print ("energy = ", a.Energy(gfu.vec))
    
        #solve linearized problem:
        a.Apply (gfu.vec, res)
        a.AssembleLinearization (gfu.vec)
        inv = a.mat.Inverse(gfu.space.FreeDofs())
        du.data = inv * res
    

        #update iteration
        gfu.vec.data -= du

        #stopping criteria
        stopcritval = sqrt(abs(InnerProduct(du,res)))
        print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)
        if stopcritval < 1e-13:
            break
        Redraw(blocking=True)
