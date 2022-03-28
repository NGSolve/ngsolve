from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import unit_cube, CSGeometry, Cylinder, Sphere, Pnt, Plane, Vec
import netgen.meshing as meshing
from ngsolve.meshes import MakeStructured2DMesh,MakeStructured3DMesh, MakeStructuredSurfaceMesh
from space_utils import *

def Test(mesh, space, order, idop=lambda cf : cf, trace=None, ttrace=None, diffops=None, vb=VOL, set_dual=[False], addorder=0, sym=False, dev=False, facet=False, **kwargs):
    fes = space(mesh, order=order+addorder, dim=kwargs.get("dim", 1))
    gf = GridFunction(fes)
    print("space = ", space)
    print("gf.dims = ", gf.dims)
    cf = GetDiffOp("id", order, dim=mesh.dim, dims=gf.dims, sym=sym, dev=dev, vb=vb)
    print("cf dims = ", cf.dims)

    dx_vol  = mesh.Materials(".*")   if vb==VOL else mesh.Boundaries(".*")
    dx_bnd  = mesh.Boundaries(".*")  if vb==VOL else mesh.BBoundaries(".*")
    dx_bbnd = mesh.BBoundaries(".*") if vb==VOL else mesh.BBBoundaries(".*")

    if ttrace:
        gf.Set(cf, dual=False, definedon=dx_bbnd)
        assert sqrt(Integrate(InnerProduct(gf-ttrace(cf),gf-ttrace(cf)), mesh, definedon=dx_bbnd)) < 1e-10
    if trace:
        gf.Set(cf, dual=False, definedon=dx_bnd)
        assert sqrt(Integrate(InnerProduct(gf-trace(cf),gf-trace(cf)), mesh, definedon=dx_bnd)) < 1e-10
    if idop:
        for dual in set_dual:
            gf.Set(cf, dual=dual, definedon=dx_vol)
            assert sqrt(Integrate(InnerProduct(gf-idop(cf),gf-idop(cf))*dx(definedon=dx_vol,element_boundary=facet), mesh)) < 1e-10
    if diffops:
        for diffop in diffops:
            cfdiffop = GetDiffOp(diffop.lower(), order, dim=mesh.dim, dims=gf.dims, sym=sym, dev=dev, vb=vb)
            try:
                gfdiffop = globals()[diffop](gf)
            except:
                gfdiffop = None
            if gfdiffop:
                assert sqrt(Integrate( InnerProduct(gfdiffop-cfdiffop,gfdiffop-cfdiffop), mesh, definedon=dx_vol)) < 5e-8
            else:
                assert sqrt(Integrate( InnerProduct(gf.Operator(diffop,vb)-cfdiffop,gf.Operator(diffop,vb)-cfdiffop), mesh, definedon=dx_vol)) < 5e-8
    return


def test_fespaces_2d():
    n_2d = specialcf.normal(2)
    Ptau_2d = Id(2) - OuterProduct(n_2d,n_2d)
    Pn_2d = OuterProduct(n_2d,n_2d)

    # unstructured trig mesh
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3,quad_dominated=False))
    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["hesse", "Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["Grad"], vb=VOL, set_dual=[True,False], dim=2)
    Test(mesh=mesh, space=VectorH1, order=2, trace = lambda cf : cf, diffops=["hesse", "div", "Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=L2, order=2, diffops=["Grad","hesse"], vb=VOL, set_dual=[False])
    Test(mesh=mesh, space=VectorL2, order=2, diffops=["Grad"], vb=VOL, set_dual=[False])
    Test(mesh=mesh, space=HCurl, order=2, trace = lambda cf : Ptau_2d*cf, diffops=["curl","Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=HDiv, order=2, trace = lambda cf : Pn_2d*cf, diffops=["div","Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=HDivDiv, order=2, trace = lambda cf : Pn_2d*cf*Pn_2d, diffops=["div"], vb=VOL, set_dual=[False], sym=True)
    Test(mesh=mesh, space=HCurlCurl, order=2, trace = lambda cf : Ptau_2d*cf*Ptau_2d, diffops=["curl","inc", "christoffel","christoffel2"], vb=VOL, set_dual=[False], sym=True)
    Test(mesh=mesh, space=HCurlDiv, order=2, trace = lambda cf : Ptau_2d*cf*Pn_2d, diffops=["div","curl"], vb=VOL, set_dual=[False], dev=True)
    Test(mesh=mesh, space=FacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True)
    Test(mesh=mesh, space=FacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True, dim=2)
    Test(mesh=mesh, space=VectorFacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True)
    Test(mesh=mesh, space=NormalFacetFESpace, order=2, trace = lambda cf : Pn_2d*cf, vb=VOL, set_dual=[], facet=True)# dual=True: netgen.libngpy._meshing.NgException: normal-facet element evaluated not at BND
    #Test(mesh=mesh, space=TangentialFacetFESpace, order=2, idop = lambda cf : Ptau_2d*cf, trace = lambda cf : Ptau_2d*cf, vb=VOL, set_dual=[True], facet=True) #idop not Ptau*cf?
    
    """
    # unstructured (= non-affine) quad mesh
    mesh = MakeStructured2DMesh(quads=True, nx=3,ny=3, mapping = lambda x,y : (1.3*x*(0.4+0.4*y)**2, 0.75*y))
    #Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["hesse", "Grad"], vb=VOL, set_dual=[True,False], addorder=1)
    
    # unstructured trig/quad mixed mesh (boundary could be curved?)
    geo = SplineGeometry()
    geo.AddRectangle((0,0), (1,1))
    geo.AddCircle( (0.7,0.5), r=0.1, leftdomain=0, rightdomain=1)
    mesh = Mesh(geo.GenerateMesh(quad_dominated=True,maxh=0.5))
    #Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["hesse", "Grad"], vb=VOL, set_dual=[True,False], addorder=1)
    """
    return

def test_fespaces_3d():
    n_3d = specialcf.normal(3)
    t_3d = specialcf.tangential(3)
    Ptau_3d = Id(3) - OuterProduct(n_3d,n_3d)
    Pn_3d = OuterProduct(n_3d,n_3d)
    Pt_3d = OuterProduct(t_3d,t_3d)

    # unstructured tet mesh
    mesh = MakeStructured3DMesh(hexes=False, nx=2,ny=2,nz=2, prism=False, mapping = lambda x,y,z : (x*(0.4+0.4*y)**2,0.75*y,1.25*z) )

    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["hesse", "Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["Grad"], vb=VOL, set_dual=[True,False], dim=3)
    Test(mesh=mesh, space=VectorH1, order=2, trace = lambda cf : cf, diffops=["hesse", "div", "Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=L2, order=2, diffops=["Grad","hesse"], vb=VOL, set_dual=[False])
    Test(mesh=mesh, space=VectorL2, order=2, diffops=["Grad"], vb=VOL, set_dual=[False])
    Test(mesh=mesh, space=HCurl, order=2, trace = lambda cf : Ptau_3d*cf, diffops=["curl","Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=HDiv, order=2, trace = lambda cf : Pn_3d*cf, diffops=["div","Grad"], vb=VOL, set_dual=[True,False])
    Test(mesh=mesh, space=HDivDiv, order=2, trace = lambda cf : Pn_3d*cf*Pn_3d, diffops=["div"], vb=VOL, set_dual=[False], sym=True)
    Test(mesh=mesh, space=HCurlCurl, order=2, trace = lambda cf : Ptau_3d*cf*Ptau_3d, diffops=["curl","inc","christoffel","christoffel2"], vb=VOL, set_dual=[False], sym=True)
    Test(mesh=mesh, space=HCurlDiv, order=2, trace = lambda cf : Ptau_3d*cf*Pn_3d, diffops=["div"], vb=VOL, set_dual=[False], dev=True)
    Test(mesh=mesh, space=FacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True)
    Test(mesh=mesh, space=FacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True, dim=3)
    Test(mesh=mesh, space=VectorFacetFESpace, order=2, trace = lambda cf : cf, vb=VOL, set_dual=[True], facet=True)
    Test(mesh=mesh, space=NormalFacetFESpace, order=2, trace = lambda cf : Pn_3d*cf, vb=VOL, set_dual=[], facet=True)
    
    
    """
    # other unstructured tet mesh
    #geo = CSGeometry()
    #geo.Add(Sphere(Pnt(0,0,0), 1))    
    #mesh = Mesh(geo.GenerateMesh(maxh=0.125))
    #Test(mesh=mesh, order=2, diffops=True, vb=VOL, set_dual=[True,False])
    
    # unstructured hex mesh
    mesh = MakeStructured3DMesh(hexes=True, nx=2,ny=2,nz=2, prism=False, mapping = lambda x,y,z : (x*(0.4+0.4*y)**2,0.75*y,1.25*z) )
    Test(mesh=mesh, order=2, vb=VOL, addorder=1)
    
    # unstructured prism mesh
    mesh = MakeStructured3DMesh(hexes=False, nx=2,ny=2,nz=2, prism=True, mapping = lambda x,y,z : (x*(0.4+0.4*y)**2,0.75*y,1.25*z) )
    Test(mesh=mesh, order=2, vb=VOL)
    """
    
    return


def test_fespaces_surface():
    n_3d = specialcf.normal(3)
    t_3d = specialcf.tangential(3)
    Ptau_3d = Id(3) - OuterProduct(n_3d,n_3d)
    Pn_3d = OuterProduct(n_3d,n_3d)
    Pt_3d = OuterProduct(t_3d,t_3d)

    # unstructured trig surface mesh (surface could be curved?)
    mesh = MakeStructuredSurfaceMesh(quads=False, nx=3, ny=3, mapping = lambda x,y,z : ( (x-0.5), (y-0.5), (x-0.5)**2*(0.7+0.2*y)-(y-0.5)**2) )
    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["Grad", "hesseboundary"], vb=BND, set_dual=[True,False])
    Test(mesh=mesh, space=H1, order=2, trace = lambda cf : cf, diffops=["Grad"], vb=BND, set_dual=[True,False],dim=3)
    Test(mesh=mesh, space=VectorH1, order=2, trace = lambda cf : cf, diffops=["Grad","div"], vb=BND, set_dual=[True,False])
    Test(mesh=mesh, space=SurfaceL2, order=2, diffops=None, vb=BND, set_dual=[False])
    Test(mesh=mesh, space=HCurl, order=2, idop = lambda cf : Ptau_3d*cf, trace = lambda cf : Pt_3d*cf, diffops=None, vb=BND, set_dual=[True,False])
    Test(mesh=mesh, space=HCurlCurl, order=2, idop = lambda cf : Ptau_3d*cf*Ptau_3d, trace = lambda cf : Pt_3d*cf*Pt_3d, diffops=None, vb=BND, set_dual=[True,False], sym=True)
    Test(mesh=mesh, space=FacetSurface, order=2, trace = lambda cf : cf, diffops=None, vb=BND, set_dual=[True], facet=True)
    Test(mesh=mesh, space=FacetSurface, order=2, trace = lambda cf : cf, diffops=None, vb=BND, set_dual=[True], facet=True, dim=3)
    Test(mesh=mesh, space=VectorFacetSurface, order=2, trace = lambda cf : cf, diffops=None, vb=BND, set_dual=[True], facet=True)
    
    #Test(mesh=mesh, space=NormalFacetSurface, order=2, diffops=None, vb=BND, set_dual=[True], facet=True) #no dual diffop
    ##Test(mesh=mesh, space=HDivSurface, order=0, idop = lambda cf : Ptau_3d*cf, trace = None, diffops=None, vb=BND, set_dual=[True], addorder=0)
    
    """
    # unstructured quad surface mesh (surface could be curved?)
    mesh = MakeStructuredSurfaceMesh(quads=True, nx=3, ny=3, mapping = lambda x,y,z : ( (x-0.5), (y-0.5), (x-0.5)**2*(0.7+0.2*y)-(y-0.5)**2) )
    Test(mesh=mesh, order=2, vb=BND, addorder=1)
    
    # unstructured trig/quad mixed surface mesh (surface could be curved?)
    geo       = CSGeometry()
    cyl       = Cylinder(Pnt(0,0,0), Pnt(1,0,0), 1)
    bot       = Plane(Pnt(0,0,0), Vec(0,0.3,-1))
    right     = Plane( Pnt(3,0,0), Vec(1,-0.1,0.2))
    left      = Plane(Pnt(0,0,0), Vec(-1,0.2,-0.03))
    additional= Plane(Pnt(0,0,0), Vec(-1,-1,-1))
    finitecyl = cyl * bot * left * right*additional
    geo.AddSurface(cyl, finitecyl)
    mesh = Mesh(geo.GenerateMesh(maxh=0.5, quad_dominated=True))
    Test(mesh=mesh, order=2, vb=BND, addorder=1)
    """
    return

    
if __name__ == "__main__":
    test_fespaces_2d()
    test_fespaces_3d()
    test_fespaces_surface()
