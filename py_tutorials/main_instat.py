from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *
from ngsolve.solve import *

from math import sin
from time import sleep


mesh = Mesh("square.vol.gz")

v = FESpace ("h1ho", mesh, order=5)  # , dirichlet=[1,2])
u = GridFunction (v, name="potential")
Draw (u, sd=1)

f = LinearForm (v)
# f.Add (LFI (name = "source", dim = 2, coef = 1))
f.Assemble()

a = BilinearForm (v, symmetric = True) 
a.Add (BFI ("laplace", dim = 2, coef = 0.01))
a.Assemble()

m = BilinearForm (v, symmetric = True)
m.Add (BFI ("mass", dim = 2, coef = 1))
m.Assemble()


tau = 0.01

mstar = a.mat.CreateMatrix()
mstar.AsVector().data = tau * a.mat.AsVector() + m.mat.AsVector()
inv = mstar.Inverse(v.FreeDofs())
        
d = u.vec.CreateVector()
w = u.vec.CreateVector()

x = CoordCF(0)
y = CoordCF(1)
u.Set (exp(-40*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))))

Redraw(blocking=True)

# print ("Switch visualization to 'Solution', and set 'ScalarFunction' to 'gfu'")
input("I'm waiting for you ... ")

t = 0.0;
while (t <= 1):
    t += tau
    d.data = sin(t) * f.vec - a.mat * u.vec
    w.data = inv * d
    u.vec.data += tau * w
    
    print ("t = ", t)
    Redraw(blocking=True)
    sleep (0.01)
           
      
