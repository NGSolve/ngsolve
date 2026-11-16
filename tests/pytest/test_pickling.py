from netgen.geom2d import *
from netgen.csg import unit_cube
from ngsolve import *
import pickle
import numpy
import io, pytest

def test_pickle_volume_fespaces():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    spaces = [H1(mesh,order=3,dirichlet=[1,2,3,4]), VectorH1(mesh,order=3,dirichlet=[1,2,3,4]), L2(mesh,order=3), VectorL2(mesh,order=3), SurfaceL2(mesh,order=3), HDivDiv(mesh,order=3,dirichlet=[1,2,3,4]), TangentialFacetFESpace(mesh,order=3,dirichlet=[1,2,3,4]), FacetFESpace(mesh,order=3,dirichlet=[1,2,3,4]), NumberSpace(mesh), HDiv(mesh,order=3,dirichlet=[1,2,3,4]), HCurl(mesh,order=3,dirichlet=[1,2,3,4])]
    
    for space in spaces:
        data = pickle.dumps(space)
        space2 = pickle.loads(data)
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

    spaces = [HDivDivSurface(mesh,order=3,dirichlet=[1,2,3,4]), FacetSurface(mesh,order=3,dirichlet=[1,2,3,4]), HDivSurface(mesh,order=3,dirichlet=[1,2,3,4])]

    for space in spaces:
        data = pickle.dumps(space)
        space2 = pickle.loads(data)
        assert space.ndof == space2.ndof
        
def test_pickle_gridfunction_real():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh,order=3,dirichlet=[1,2,3,4])
    u = GridFunction(fes,"u")
    u.Set(x*y)

    data = pickle.dumps((u, grad(u)))
    u2, gradu2 = pickle.loads(data)
    assert sqrt(Integrate((u-u2)*(u-u2),mesh)) < 1e-14
    assert sqrt(Integrate(InnerProduct(grad(u)-gradu2, grad(u)-gradu2), mesh)) < 1e-14

def test_pickle_multidim():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.4))
    fes = H1(mesh,order=3,dim=3)
    u = GridFunction(fes)
    u.Set((1,2,3))
    pickled = pickle.dumps((u,grad(u)))
    u2, gradu2 = pickle.loads(pickled)
    assert numpy.linalg.norm(u.vec.FV().NumPy() - u2.vec.FV().NumPy()) < 1e-14
    assert Integrate(Norm(grad(u)-gradu2), mesh) < 1e-14

def test_pickle_secondorder_mesh():
    m = unit_cube.GenerateMesh(maxh=0.4)
    m.SecondOrder()
    mesh = Mesh(m)
    fes = H1(mesh, order=3)
    u = GridFunction(fes)
    u.Set(x)
    pickled = pickle.dumps(u)
    u2 = pickle.loads(pickled)
    assert numpy.linalg.norm(u.vec.FV().NumPy() - u2.vec.FV().NumPy()) < 1e-14


def test_pickle_gridfunction_complex():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = HCurl(mesh,order=3,complex=True)
    u = GridFunction(fes,"u")
    u.Set((x*y,1J * y))
    data = pickle.dumps((u, curl(u)))
    u2, curlu2  = pickle.loads(data)
    difvec = u.vec.CreateVector()
    difvec.data = u.vec - u2.vec
    assert Norm(difvec) < 1e-12
    assert Integrate(Norm(curl(u)-curlu2), mesh) < 1e-14

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

def test_pickle_allsamespace():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))

    fes1 = H1(mesh, order=1)

    fes = fes1**3

    gfu = GridFunction(fes)
    gfu.Set ( (x,y,z) )
    data = pickle.dumps(gfu)
    gfu2 = pickle.loads(data)
    assert Integrate(Norm(gfu - gfu2), mesh) < 1e-14

@pytest.mark.parametrize("phase", [None, [-1]*10, [1J]*10])
def test_pickle_periodic(phase):
    periodic = SplineGeometry()
    pnts = [ (0,0), (1,0), (1,1), (0,1) ]
    pnums = [periodic.AppendPoint(*p) for p in pnts]
    periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
    lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
    periodic.Append ( ["line", pnums[2], pnums[3]], bc="outer")
    periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
    mesh = Mesh(periodic.GenerateMesh(maxh=0.2))
    is_complex = (phase is not None) and (type(phase[0]) != float)
    fes = Periodic(H1(mesh,order=3, complex=is_complex, dirichlet="outer"),phase=phase)
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
    data = pickle.dumps(u)
    u2 = pickle.loads(data)
    assert u2.space.is_complex == is_complex
    assert Integrate(Norm(u-u2),mesh) < 1e-14


def test_pickle_compress():
    periodic = SplineGeometry()
    pnts = [ (0,0), (1,0), (1,1), (0,1) ]
    pnums = [periodic.AppendPoint(*p) for p in pnts]
    periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
    lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
    periodic.Append ( ["line", pnums[2], pnums[3]], bc="outer")
    periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
    mesh = Mesh(periodic.GenerateMesh(maxh=0.2))
    fes = Compress(Periodic(H1(mesh,order=3, dirichlet="outer")))
    data = pickle.dumps(fes)
    fes2 = pickle.loads(data)
    assert fes.ndof == fes2.ndof

def test_pickle_CoefficientFunctions():
    import netgen.csg as csg
    import numpy as np
    import numpy.testing as npt
    mesh = Mesh(csg.unit_cube.GenerateMesh(maxh=1000))
    cfs = []
    cfs.append(CoefficientFunction(1)) # ConstantCF
    cfs.append(CoefficientFunction(1J)) # ConstantCFC
    cfs.append(Parameter(3)) # ParameterCF
    cfs.append(CoefficientFunction([cfs[2]])) # Domainwise with parameter inside
    # coordinate
    cfs.append(x)
    cfs.append(y)
    cfs.append(z)
    # unary
    cfs.append(sin(x))
    cfs.append(cos(y))
    cfs.append(floor(z))
    # binary
    cfs.append(x**y)
    cfs.append(y*z)
    cfs.append(x/y)
    cfs.append(x+cfs[2])
    cfs.append(cfs[2]-y)
    vec = CoefficientFunction((cfs[2],1,2))
    mat = CoefficientFunction((5,cfs[2],0,2),dims=(2,2))
    cfs.append(vec) # Vectorial
    cfs.append(vec[0]) # Component
    cfs.append(2 * vec) # scale
    cfs.append(vec * vec) # MultVecVec
    cfs.append(Norm(vec)) # Norm
    cfs.append(mat) # Mat
    cfs.append(mat.trans)
    cfs.append(mat.Eig())
    cfs.append(IfPos(x-cfs[2], y, z)) # Ifpos
    imagcf = x * 1J + y
    cfs.append(imagcf)
    cfs.append(imagcf.real)
    cfs.append(imagcf.imag)
    cfs.append(Norm(imagcf))
    to_compile = Inv(Cof(mat))
    cfs.append(to_compile.Compile())
    cfs.append(to_compile.Compile(realcompile=True, wait=True))
    cfs.append(to_compile.Compile(realcompile=True, wait=False))
    mat55 = CoefficientFunction((5,cfs[2],0,2, 1,
                                 0,     1,0,0,0,
                                 0,     0,2,1,0,
                                 1,     0,0,1,0,
                                 0,     2,0,4,1),
                                dims=(5,5))
    cfs.append(Inv(mat55))
    cfs.append(mat55[2,1:])
    def _DetES(A, **options):
        return fem.Einsum('ijk,i,j,k->', fem.LeviCivitaSymbol(3), *[A[i, :] for i in range(3)], **options)
    mat33 = mat55[:3,:3]
    cfs.append(_DetES(mat33))
    cfs.append(_DetES(mat33, optimize_path=True))

    mp = mesh(0.3,0.4,0.5)
    compiled_vals = to_compile(mp)
    dump = pickle.dumps(cfs)
    lcfs = pickle.loads(dump)
    assert lcfs[0](mp) == 1
    assert lcfs[1](mp) == 1J
    assert lcfs[2].Get() == 3
    assert lcfs[3](mp) == 3, lcfs[3](mp)
    lcfs[2].Set(5)
    assert lcfs[2](mp) == 5
    assert lcfs[3](mp) == 5
    npt.assert_approx_equal(lcfs[4](mp),0.3)
    npt.assert_approx_equal(lcfs[5](mp),0.4)
    npt.assert_approx_equal(lcfs[6](mp),0.5)
    npt.assert_approx_equal(lcfs[7](mp),sin(0.3))
    npt.assert_approx_equal(lcfs[8](mp),cos(0.4))
    assert lcfs[9](mp) == 0
    npt.assert_approx_equal(lcfs[10](mp), 0.3**0.4)
    npt.assert_approx_equal(lcfs[11](mp), 0.4*0.5)
    npt.assert_approx_equal(lcfs[12](mp), 0.3/0.4)
    npt.assert_approx_equal(lcfs[13](mp), 5.3)
    npt.assert_approx_equal(lcfs[14](mp), 4.6)
    assert lcfs[15](mp) == (5,1,2)
    lcfs[2].Set(8)
    assert lcfs[16](mp) == 8
    assert lcfs[17](mp) == (16,2,4)
    assert lcfs[18](mp) == 69
    assert lcfs[19](mp) == sqrt(69)
    assert lcfs[20](mp) == (5,8,0,2), lcfs[20](mp)
    assert lcfs[21](mp) == (5,0,8,2)
    lcfs[2].Set(0)
    assert lcfs[22](mp) == (1,0,0,1,5,2)
    lcfs[2].Set(0.5)
    npt.assert_approx_equal(lcfs[23](mp),0.5)
    npt.assert_approx_equal(lcfs[23](mesh(0.7,0.1,0.5)),0.1)
    npt.assert_approx_equal(lcfs[24](mp).real,0.4)
    npt.assert_approx_equal(lcfs[24](mp).imag,0.3)
    npt.assert_approx_equal(lcfs[25](mp),0.4)
    npt.assert_approx_equal(lcfs[26](mp),0.3)
    npt.assert_approx_equal(lcfs[27](mp),sqrt(0.3*0.3+0.4*0.4))
    lcfs[2].Set(cfs[2].Get())
    npt.assert_array_almost_equal(lcfs[28](mp), compiled_vals)
    npt.assert_array_almost_equal(lcfs[29](mp), compiled_vals)
    npt.assert_array_almost_equal(lcfs[30](mp), compiled_vals)
    npt.assert_array_almost_equal(lcfs[31](mp), cfs[31](mp))
    npt.assert_array_almost_equal(lcfs[32](mp), cfs[32](mp))
    npt.assert_array_almost_equal(lcfs[33](mp), cfs[33](mp))
    npt.assert_array_almost_equal(lcfs[34](mp), cfs[33](mp))

def test_pickle_sparsecholesky():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    fes = H1(mesh, order=2, dirichlet=".*")
    u,v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    f = LinearForm(v * dx).Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    u = GridFunction(fes)
    u2 = GridFunction(fes)
    data = pickle.dumps(inv)
    inv2 = pickle.loads(data)
    u.vec.data = inv * f.vec
    u2.vec.data = inv2 * f.vec - u.vec
    assert Norm(u2.vec) < 1e-16


# region-wise CFs pickled by NGSolve 6.2.2604, when the class was DomainWiseCoefficientFunction
_old_domainwise_pickles = {
  "list": (
    "gAWVngMAAAAAAACMC25nc29sdmUuZmVtlIwTQ29lZmZpY2llbnRGdW5jdGlvbpSTlCmBlF2UKGgCKYGUXZQoaAIpgZRdlChoAimB"
    "lF2UKGgCKYGUXZQoQzL9////IgAAAG5nZmVtOjpDb25zdGFudENvZWZmaWNpZW50RnVuY3Rpb24AAAAAAAAUQJRDUQIAAAAAAAAA"
    "BgAAAG5ldGdlbhUAAAB2Ni4yLjI2MDQtOC1nZjNlZDM1OTAHAAAAbmdzb2x2ZRcAAAB2Ni4yLjI2MDQtMTgtZzMxNWVhMzQ4M5RD"
    "CAAAAAAAAAAAlGWFlGJDSf3///8rAAAAbmdmZW06OmNsX1VuYXJ5T3BDRjxuZ2ZlbTo6R2VuZXJpY0lkZW50aXR5PgEAAAAgAQAA"
    "AAAAAAAAAAAAAAAAAACUQ1ECAAAAAAAAAAYAAABuZXRnZW4VAAAAdjYuMi4yNjA0LTgtZ2YzZWQzNTkwBwAAAG5nc29sdmUXAAAA"
    "djYuMi4yNjA0LTE4LWczMTVlYTM0ODOUQwgAAAAAAAAAAJRlhZRiQ0n9////KwAAAG5nZmVtOjpjbF9VbmFyeU9wQ0Y8bmdmZW06"
    "OkdlbmVyaWNJZGVudGl0eT4BAAAAIAEAAAAAAAAAAAAAAAAAAAAAlENRAgAAAAAAAAAGAAAAbmV0Z2VuFQAAAHY2LjIuMjYwNC04"
    "LWdmM2VkMzU5MAcAAABuZ3NvbHZlFwAAAHY2LjIuMjYwNC0xOC1nMzE1ZWEzNDgzlEMIAAAAAAAAAACUZYWUYkNG/f///yQAAABu"
    "Z2ZlbTo6RG9tYWluV2lzZUNvZWZmaWNpZW50RnVuY3Rpb24BAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAJRDUQIAAAAAAAAABgAA"
    "AG5ldGdlbhUAAAB2Ni4yLjI2MDQtOC1nZjNlZDM1OTAHAAAAbmdzb2x2ZRcAAAB2Ni4yLjI2MDQtMTgtZzMxNWVhMzQ4M5RDCAAA"
    "AAAAAAAAlGWFlGJDSf3///8rAAAAbmdmZW06OmNsX1VuYXJ5T3BDRjxuZ2ZlbTo6R2VuZXJpY0lkZW50aXR5PgEAAAAgAQAAAAAA"
    "AAAAAAAAAAAAAACUQ1ECAAAAAAAAAAYAAABuZXRnZW4VAAAAdjYuMi4yNjA0LTgtZ2YzZWQzNTkwBwAAAG5nc29sdmUXAAAAdjYu"
    "Mi4yNjA0LTE4LWczMTVlYTM0ODOUQwgAAAAAAAAAAJRlhZRiLg=="),
  "material": (
    "gAWV6wIAAAAAAACMC25nc29sdmUuZmVtlIwTQ29lZmZpY2llbnRGdW5jdGlvbpSTlCmBlF2UKGgCKYGUXZQoaAIpgZRdlChoAimB"
    "lF2UKEMy/f///yIAAABuZ2ZlbTo6Q29uc3RhbnRDb2VmZmljaWVudEZ1bmN0aW9uAAAAAAAACECUQ1ECAAAAAAAAAAYAAABuZXRn"
    "ZW4VAAAAdjYuMi4yNjA0LTgtZ2YzZWQzNTkwBwAAAG5nc29sdmUXAAAAdjYuMi4yNjA0LTE4LWczMTVlYTM0ODOUQwgAAAAAAAAA"
    "AJRlhZRiQ0n9////KwAAAG5nZmVtOjpjbF9VbmFyeU9wQ0Y8bmdmZW06OkdlbmVyaWNJZGVudGl0eT4BAAAAIAEAAAAAAAAAAAAA"
    "AAAAAAAAlENRAgAAAAAAAAAGAAAAbmV0Z2VuFQAAAHY2LjIuMjYwNC04LWdmM2VkMzU5MAcAAABuZ3NvbHZlFwAAAHY2LjIuMjYw"
    "NC0xOC1nMzE1ZWEzNDgzlEMIAAAAAAAAAACUZYWUYkNJ/f///ysAAABuZ2ZlbTo6Y2xfVW5hcnlPcENGPG5nZmVtOjpHZW5lcmlj"
    "SWRlbnRpdHk+AQAAACABAAAAAAAAAAAAAAAAAAAAAJRDUQIAAAAAAAAABgAAAG5ldGdlbhUAAAB2Ni4yLjI2MDQtOC1nZjNlZDM1"
    "OTAHAAAAbmdzb2x2ZRcAAAB2Ni4yLjI2MDQtMTgtZzMxNWVhMzQ4M5RDCAAAAAAAAAAAlGWFlGJDSv3///8kAAAAbmdmZW06OkRv"
    "bWFpbldpc2VDb2VmZmljaWVudEZ1bmN0aW9uAQAAAAAAAAABAAAAAAEAAAAAAAAAAAAAAAAAAAAAlENRAgAAAAAAAAAGAAAAbmV0"
    "Z2VuFQAAAHY2LjIuMjYwNC04LWdmM2VkMzU5MAcAAABuZ3NvbHZlFwAAAHY2LjIuMjYwNC0xOC1nMzE1ZWEzNDgzlEMIAAAAAAAA"
    "AACUZYWUYi4="),
  "boundary": (
    "gAWVlAMAAAAAAACMC25nc29sdmUuZmVtlIwTQ29lZmZpY2llbnRGdW5jdGlvbpSTlCmBlF2UKGgCKYGUXZQoQzf9////HgAAAG5n"
    "ZmVtOjpaZXJvQ29lZmZpY2llbnRGdW5jdGlvbgEAAAAAAAAAAAAAAAAAAAAAlENRAgAAAAAAAAAGAAAAbmV0Z2VuFQAAAHY2LjIu"
    "MjYwNC04LWdmM2VkMzU5MAcAAABuZ3NvbHZlFwAAAHY2LjIuMjYwNC0xOC1nMzE1ZWEzNDgzlEMIAAAAAAAAAACUZYWUYmgFaAVo"
    "AimBlF2UKGgCKYGUXZQoaAIpgZRdlChDMv3///8iAAAAbmdmZW06OkNvbnN0YW50Q29lZmZpY2llbnRGdW5jdGlvbgAAAAAAAABA"
    "lENRAgAAAAAAAAAGAAAAbmV0Z2VuFQAAAHY2LjIuMjYwNC04LWdmM2VkMzU5MAcAAABuZ3NvbHZlFwAAAHY2LjIuMjYwNC0xOC1n"
    "MzE1ZWEzNDgzlEMIAAAAAAAAAACUZYWUYkNJ/f///ysAAABuZ2ZlbTo6Y2xfVW5hcnlPcENGPG5nZmVtOjpHZW5lcmljSWRlbnRp"
    "dHk+AQAAACABAAAAAAAAAAAAAAAAAAAAAJRDUQIAAAAAAAAABgAAAG5ldGdlbhUAAAB2Ni4yLjI2MDQtOC1nZjNlZDM1OTAHAAAA"
    "bmdzb2x2ZRcAAAB2Ni4yLjI2MDQtMTgtZzMxNWVhMzQ4M5RDCAAAAAAAAAAAlGWFlGJDSf3///8rAAAAbmdmZW06OmNsX1VuYXJ5"
    "T3BDRjxuZ2ZlbTo6R2VuZXJpY0lkZW50aXR5PgEAAAAgAQAAAAAAAAAAAAAAAAAAAACUQ1ECAAAAAAAAAAYAAABuZXRnZW4VAAAA"
    "djYuMi4yNjA0LTgtZ2YzZWQzNTkwBwAAAG5nc29sdmUXAAAAdjYuMi4yNjA0LTE4LWczMTVlYTM0ODOUQwgAAAAAAAAAAJRlhZRi"
    "Q0r9////JAAAAG5nZmVtOjpEb21haW5XaXNlQ29lZmZpY2llbnRGdW5jdGlvbgQAAAAAAAAAAQEAAAABAAAAAAAAAAAAAAAAAAAA"
    "AJRDUQIAAAAAAAAABgAAAG5ldGdlbhUAAAB2Ni4yLjI2MDQtOC1nZjNlZDM1OTAHAAAAbmdzb2x2ZRcAAAB2Ni4yLjI2MDQtMTgt"
    "ZzMxNWVhMzQ4M5RDCAAAAAAAAAAAlGWFlGIu"),
}

def test_pickle_old_domainwise_cf():
    import base64
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
    cfs = {k: pickle.loads(base64.b64decode(v)) for k, v in _old_domainwise_pickles.items()}
    assert Integrate(cfs["list"], mesh, VOL) == pytest.approx(5)
    assert Integrate(cfs["material"], mesh, VOL) == pytest.approx(3)
    assert Integrate(cfs["boundary"], mesh, BND) == pytest.approx(2)
    with pytest.raises(Exception, match="evaluated on BND"):
        Integrate(cfs["material"], mesh, BND)
    # loaded objects pickle again
    assert Integrate(pickle.loads(pickle.dumps(cfs["material"])), mesh, VOL) == pytest.approx(3)

if __name__ == "__main__":
    test_pickle_volume_fespaces()
    test_pickle_surface_fespaces()
    test_pickle_gridfunction_real()
    test_pickle_gridfunction_complex()
    test_pickle_compoundfespace()
    test_pickle_hcurl()
    test_pickle_periodic()
    test_pickle_CoefficientFunctions()
    test_pickle_multidim()
    test_pickle_secondorder_mesh()
    test_pickle_compress()