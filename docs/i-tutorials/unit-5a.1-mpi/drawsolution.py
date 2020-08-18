from ngsolve import *

mesh = Mesh("square.vol")
fes = H1(mesh, order=3)
gfu = GridFunction(fes)
gfu.Load("solution.sol", parallel=True)

Draw (gfu)
