from ngsolve import *
from ngsolve.meshes import Make1DMesh
ngsmesh = Make1DMesh(10)

V = H1(ngsmesh, order=2, dirichlet="left")
print ("freedofs:\n", V.FreeDofs())


a = BilinearForm(V)
a += Laplace(1)
a.Assemble()

print ("mat = \n", a.mat)

f = LinearForm(V)    
f += Source(1)
# f += Neumann(1)
f.Assemble()

print ("rhs = \n", f.vec)

u = GridFunction(V)
u.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec

print ("sol =\n", u.vec)


print ("u(0.5) =", u(0.5))



pnts = []
for i in range(101): pnts.append (i/100)

pnts_vals = [ (x,u(x)) for x in pnts if ngsmesh.Contains(x)]
# print (pnts_vals)


import matplotlib.pyplot as plt
pnts,vals = zip(*pnts_vals)
plt.plot(pnts,vals, "-*")

plt.ion()
plt.show()


