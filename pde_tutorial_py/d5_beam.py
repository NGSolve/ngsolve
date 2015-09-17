from ngsolve import *

mesh = Mesh("beam.vol")

V = H1(mesh, order=3, dim=3, dirichlet=[1])
u = GridFunction(V)

E = 2.1E11
nu = 0.2
a = BilinearForm(V)
a += BFI("elasticity", coef=[E,nu])

c = Preconditioner(a,"multigrid")
# c = Preconditioner(a,"direct")
a.Assemble()


traction = DomainConstantCF ([0,1e5])
f = LinearForm(V)
f += BlockLFI (Neumann(traction), dim=3, comp=2)
f.Assemble()


# u.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec

BVP(bf=a, gf=u, lf=f, pre=c, maxsteps=1000).Do()

Draw (u)


