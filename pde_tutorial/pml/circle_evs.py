from ngsolve import *
from netgen.geom2d import *
from matplotlib.pyplot import plot,show
from numpy import sqrt


def ccc(r1,r2,r3):
    geom = SplineGeometry()
    geom.AddCircle((0,0),r1,leftdomain=0,rightdomain=1)
    geom.AddCircle((0,0),r2,leftdomain=1,rightdomain=2)
    geom.AddCircle((0,0),r3,leftdomain=2,rightdomain=0)
    return geom

def cbb(r1,r2,r3):
    geom = SplineGeometry()
    geom.AddCircle((0,0),r1,leftdomain=0,rightdomain=1)
    geom.AddRectangle((-r2,-r2),(r2,r2),leftdomain=1,rightdomain=2)
    geom.AddRectangle((-r3,-r3),(r3,r3),leftdomain=2,rightdomain=0)
    return geom

def circle_evs(geom,pmltrafo):
    # mesh = Mesh("cavity.vol.gz")
    mesh = Mesh(geom.GenerateMesh (maxh=0.05))
    print ("nv = ", mesh.nv)
    mesh.Curve(5)
    # define constant hpref = 4

    mesh.SetPML(pmltrafo,1)
    #SetPMLParameters (rad=0.8, alpha=2)

    fes = H1(mesh, order=4, complex=True, dirichlet=[3])
    u = fes.TrialFunction()
    v = fes.TestFunction()

    a = BilinearForm (fes, symmetric=True)
    a += SymbolicBFI(grad(u)*grad(v))

    m = BilinearForm (fes, symmetric=True)
    m += SymbolicBFI(u*v)
    a.Assemble()
    m.Assemble()

    u = GridFunction(fes, multidim=50)
    lam = ArnoldiSolver(a.mat, m.mat, fes.FreeDofs(), u.vecs, shift=1)
    print ("lam: ", lam)
    Draw(u)
    return lam
geom = ccc(0.5,0.6,1.2)
pr=pml.Radial(0.6,2j)
lamr=circle_evs(geom,pr)

geom = cbb(0.5,0.6,1.2)
pc=pml.Cartesian((-0.6,-0.6),(0.6,0.6),2j)
lamc=circle_evs(geom,pc)


plot([sqrt(l).real for l in lamr],[sqrt(l).imag for l in lamr],'xr')
plot([sqrt(l).real for l in lamc],[sqrt(l).imag for l in lamc],'xb')
show()

