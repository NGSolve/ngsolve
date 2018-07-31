
from netgen.csg import *
from ngsolve import *
from numpy import pi


def test_periodicH1():
    # incoming wave direction
    d=[1./sqrt(3),1./sqrt(3),1/sqrt(3)]
    
    # wave number
    k = 10
    Lx = 2*pi/(k*d[0])
    Ly = 2*pi/(k*d[1])
    Lz = 2*pi/(k*d[2])
    print(Lx)

    geo = CSGeometry()

    left = Plane(Pnt(0,0,0),Vec(-1,0,0))
    right = Plane(Pnt(Lx,0,0),Vec(1,0,0))
    back = Plane(Pnt(0,0,0),Vec(0,-1,0))
    front = Plane(Pnt(0,Ly,0),Vec(0,1,0))
    bot = Plane(Pnt(0,0,0),Vec(0,0,-1))
    top = Plane(Pnt(0,0,Lz),Vec(0,0,1))

    
    cube = left * right * top * bot * back * front
    
    geo.Add(cube)
    geo.PeriodicSurfaces(left,right)
    geo.PeriodicSurfaces(back,front)
    geo.PeriodicSurfaces(bot,top)

    mesh = Mesh(geo.GenerateMesh(maxh=min(Lx,Ly,Lz)/10))
    
    fes = Periodic(H1(mesh,order=2, complex=True))
    
    u,v = fes.TrialFunction(), fes.TestFunction()
    
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v) + u*v)

    c = Preconditioner(a,"direct")

    f = LinearForm(fes)
    f += SymbolicLFI((k*k+1) * exp(1J * k * (d[0] * x + d[1] * y + d[2] * z)) * v)

    u = GridFunction(fes)

    u_exact = exp(1J*k*(d[0] * x + d[1] * y + d[2] * z))

    with TaskManager():
        a.Assemble()
        f.Assemble()
        ainv = CGSolver(a.mat,c.mat)
        u.vec.data = ainv * f.vec
        error = sqrt(Integrate(Conj(u-u_exact)*(u-u_exact),mesh).real)
        assert error < 1e-2

def test_quasiperiodic():
    k=1.
    Lx=1.
    Ly=2.
    Lz=1.5
    d=[1./sqrt(3),1./sqrt(3),1./sqrt(3)]

    geo = CSGeometry()

    left = Plane(Pnt(0,0,0),Vec(-1,0,0))
    right = Plane(Pnt(Lx,0,0),Vec(1,0,0))
    bot = Plane(Pnt(0,0,0),Vec(0,0,-1))
    top = Plane(Pnt(0,0,Lz),Vec(0,0,1))
    back = Plane(Pnt(0,0,0),Vec(0,-1,0))
    front = Plane(Pnt(0,Ly,0),Vec(0,1,0))

    cube = left * right * top * bot * back * front

    geo.Add(cube)
    geo.PeriodicSurfaces(left,right)
    geo.PeriodicSurfaces(back,front)
    geo.PeriodicSurfaces(bot,top)

    mesh = Mesh(geo.GenerateMesh(maxh=0.3))
    fx = exp(1J*k*d[0]*Lx)
    fy = exp(1J*k*d[1]*Ly)
    fz = exp(1J*k*d[2]*Lz)
    factors = [fx,fy,fz]
    fes = Periodic(H1(mesh,order=5,complex=True),phase=factors)
    u,v = fes.TnT()

    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v) + u * v)
    f = LinearForm(fes)
    f += SymbolicLFI((k**2+1)*exp(1J*k*(d[0]*x+d[1]*y+d[2]*z)) * v)

    u = GridFunction(fes,name="u")

    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    u_true = exp(1J * k * (d[0] * x + d[1] * y + d[2] * z))
    assert sqrt(Integrate(Conj(u_true-u)*(u_true-u),mesh).real) < 1e-8
