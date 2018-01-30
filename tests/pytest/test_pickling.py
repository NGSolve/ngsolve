from netgen.geom2d import *
from ngsolve import *
import pickle
import io

#HDivDivSurface, FacetSurface, HDivSurface
def test_pickle_volume_fespaces():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    spaces = [H1(mesh,order=3,dirichlet=[1,2,3,4]),VectorH1(mesh,order=3,dirichlet=[1,2,3,4]),L2(mesh,order=3),VectorL2(mesh,order=3),SurfaceL2(mesh,order=3), HDivDiv(mesh,order=3,dirichlet=[1,2,3,4]),VectorFacet(mesh,order=3,dirichlet=[1,2,3,4]),FacetFESpace(mesh,order=3,dirichlet=[1,2,3,4]),HybridDGSpace(mesh,order=3,dirichlet=[1,2,3,4]),NumberSpace(mesh),HDiv(mesh,order=3,dirichlet=[1,2,3,4]),HCurl(mesh,order=3,dirichlet=[1,2,3,4])]
    
    for space in spaces:
        with io.BytesIO() as f:
            pickler = pickle.Pickler(f)
            pickler.dump(space)
            data = f.getvalue()

        with io.BytesIO(data) as f:
            unpickler = pickle.Unpickler(f)
            space2 = unpickler.load()

        assert space.ndof == space2.ndof

        
def test_pickle_surface_fespaces():
    import netgen.meshing as meshing
    import netgen.csg as csg 
    geo = csg.CSGeometry()
    bottom   = csg.Plane (csg.Pnt(0,0,0), csg.Vec(0,0, 1) )
    surface = csg.SplineSurface(bottom)
    pts = [(0,0,0),(0,1,0),(1,1,0),(1,0,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for p1,p2,bc in [(0,1,"left"), (1, 2,"top"),(2,3,"right"),(3,0,"bottom")]:
        surface.AddSegment(geopts[p1],geopts[p2],bc)
    geo.AddSplineSurface(surface)
    mesh = Mesh(geo.GenerateMesh(maxh=0.3, perfstepsend=meshing.MeshingStep.MESHSURFACE))

    spaces = [HDivDivSurface(mesh,order=3,dirichlet=[1,2,3,4]),FacetSurface(mesh,order=3,dirichlet=[1,2,3,4]),HDivSurface(mesh,order=3,dirichlet=[1,2,3,4])]

    for space in spaces:
        with io.BytesIO() as f:
            pickler = pickle.Pickler(f)
            pickler.dump(space)
            data = f.getvalue()

        with io.BytesIO(data) as f:
            unpickler = pickle.Unpickler(f)
            space2 = unpickler.load()

        assert space.ndof == space2.ndof
        
def test_pickle_gridfunction_real():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh,order=3,dirichlet=[1,2,3,4])
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v))

    f = LinearForm(fes)
    f += SymbolicLFI(1*v)

    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    with io.BytesIO() as f:
        pickler = pickle.Pickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = pickle.Unpickler(f)
        u2 = unpickler.load()

    assert sqrt(Integrate((u-u2)*(u-u2),mesh)) < 1e-14


def test_pickle_gridfunction_complex():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh,order=3,complex=True,dirichlet=[1,2,3,4])
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u) * grad(v))

    f = LinearForm(fes)
    f += SymbolicLFI(1j*v)

    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    with io.BytesIO() as f:
        pickler = pickle.Pickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = pickle.Unpickler(f)
        u2 = unpickler.load()
    error = sqrt(Integrate(Conj(u-u2)*(u-u2),mesh))
    assert error.real < 1e-14 and error.imag < 1e-14

def test_pickle_hcurl():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = HCurl(mesh,order=3)
    u = GridFunction(fes)
    u.Set(CoefficientFunction((1,1)))
    with io.BytesIO() as f:
        pickler = pickle.Pickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = pickle.Unpickler(f)
        u2 = unpickler.load()
    error = sqrt(Integrate(Conj(u-u2)*(u-u2),mesh))
    assert error.real < 1e-14 and error.imag < 1e-14

def test_pickle_compoundfespace():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes1 = HDiv(mesh,order=2)
    fes2 = L2(mesh,order=1)
    fes = FESpace([fes1,fes2])
    sigma, u = fes.TrialFunction()
    tau,v = fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(sigma * tau + div(sigma)*v + div(tau)*u - 1e-10*u*v)

    f = LinearForm(fes)
    f += SymbolicLFI(-v)

    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    with io.BytesIO() as f:
        pickler = pickle.Pickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = pickle.Unpickler(f)
        u2 = unpickler.load()
    flux1, flux2 = u.components[0], u2.components[0]
    sol1, sol2 = u.components[1], u2.components[1]
    error = sqrt(Integrate((sol1-sol2)*(sol1-sol2),mesh))
    errorflux = sqrt(Integrate((flux1[0]-flux2[0])*(flux1[0]-flux2[0])+(flux1[1]-flux2[1])*(flux1[1]-flux2[1]),mesh))
    assert error < 1e-14 and errorflux < 1e-14

def test_pickle_periodic():
    periodic = SplineGeometry()
    pnts = [ (0,0), (1,0), (1,1), (0,1) ]
    pnums = [periodic.AppendPoint(*p) for p in pnts]
    periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
    lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
    periodic.Append ( ["line", pnums[2], pnums[3]], bc="outer")
    periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
    mesh = Mesh(periodic.GenerateMesh(maxh=0.2))
    fes = Periodic(H1(mesh,order=3,dirichlet="outer"))
    u,v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes,symmetric=True)
    a += SymbolicBFI(grad(u) * grad(v))
    f = LinearForm(fes)
    f += SymbolicLFI(x*v)
    u = GridFunction(fes,"u")
    with TaskManager():
        a.Assemble()
        f.Assemble()
        u.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    with io.BytesIO() as f:
        pickler = pickle.Pickler(f)
        pickler.dump(u)
        data = f.getvalue()

    with io.BytesIO(data) as f:
        unpickler = pickle.Unpickler(f)
        u2 = unpickler.load()

    assert sqrt(Integrate((u-u2)*(u-u2),mesh)) < 1e-14


if __name__ == "__main__":
    test_pickle_volume_fespaces()
    test_pickle_surface_fespaces()
    test_pickle_gridfunction_real()
    test_pickle_gridfunction_complex()
    test_pickle_compoundfespace()
    test_pickle_hcurl()
    test_pickle_periodic()
