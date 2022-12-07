from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import unit_cube, CSGeometry, Cylinder, Sphere, Pnt, Plane, Vec
import netgen.meshing as meshing
from ngsolve.meshes import MakeStructured2DMesh,MakeStructured3DMesh, MakeStructuredSurfaceMesh
from space_utils import *

def Test(mesh, space, order, idop=lambda cf : cf, trace=None, ttrace=None, diffops=None, vb=VOL, set_dual=[False], addorder=0, sym=False, dev=False, facet=False, **kwargs):
    fes = space(mesh, order=order+addorder, dim=kwargs.get("dim", 1))
    gf = GridFunction(fes)

    cf = GetDiffOp("id", order, dim=mesh.dim, dims=gf.dims, sym=sym, dev=dev, vb=vb)
    
    dx_vol  = mesh.Materials(".*")   if vb==VOL else mesh.Boundaries(".*")
    dx_bnd  = mesh.Boundaries(".*")  if vb==VOL else mesh.BBoundaries(".*")
    dx_bbnd = mesh.BBoundaries(".*") if vb==VOL else mesh.BBBoundaries(".*")

    if ttrace:
        gf.Set(cf, dual=False, definedon=dx_bbnd)
        assert sqrt(Integrate(InnerProduct(gf-ttrace(cf),gf-ttrace(cf)), mesh, definedon=dx_bbnd)) < 1e-10
    if trace:
        gf.Set(cf, dual=False, definedon=dx_bnd)
        assert sqrt(Integrate(InnerProduct(idop(gf)-trace(cf),idop(gf)-trace(cf)), mesh, definedon=dx_bnd)) < 1e-10
    if idop:
        for dual in set_dual:
            gf.Set(cf, dual=dual, definedon=dx_vol)
            assert sqrt(Integrate(InnerProduct(idop(gf)-idop(cf),idop(gf)-idop(cf))*dx(definedon=dx_vol,element_boundary=facet), mesh)) < 1e-10
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



def TestCurvatureDiffOps2D(mesh, order):
    g = CF( (1 + (x - 1/3*x**3)**2, (x - 1/3*x**3)*(y - 1/3*y**3), (x - 1/3*x**3)*(y - 1/3*y**3), 1 + (y - 1/3*y**3)**2), dims=(2,2) )
    Gauss_ex = 81*(-x**2 + 1)*(-y**2 + 1)/(9 + x**2*(x**2 - 3)**2 + y**2*(y**2 - 3)**2)**2
    Scalar_ex = 2*Gauss_ex

    Ricci_ex = -CF( ((-1 - x**2 + 2/3*x**4 - 1/9*x**6)*(x**2 - 1)*(y**2 - 1)*1/(1 + 1/9*x**6 + 1/9*y**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2)**2, (1 - 1/3*y**2)*y*(-1 + 1/3*x**2)*x*(x**2 - 1)*(y**2 - 1)*1/(1 + 1/9*x**6 + 1/9*y**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2)**2, (1 - 1/3*y**2)*y*(-1 + 1/3*x**2)*x*(x**2 - 1)*(y**2 - 1)*1/(1 + 1/9*x**6 + 1/9*y**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2)**2,(-1 - y**2 + 2/3*y**4 - 1/9*y**6)*(x**2 - 1)*(y**2 - 1)*1/(1 + 1/9*x**6 + 1/9*y**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2)**2),dims=(2,2) )
    
    fes = HCurlCurl(mesh,order=order)
    gf = GridFunction(fes)
    gf.Set(g)
    
    assert sqrt(Integrate( InnerProduct(gf.Operator("scalar")-Scalar_ex,gf.Operator("scalar")-Scalar_ex), mesh)) < 1e-8
    assert sqrt(Integrate( InnerProduct(gf.Operator("Ricci")-Ricci_ex,gf.Operator("Ricci")-Ricci_ex), mesh)) < 1e-8
    assert sqrt(Integrate( InnerProduct(gf.Operator("Einstein"),gf.Operator("Einstein")), mesh)) < 1e-8
    return

def TestCurvatureDiffOps3D(mesh, order):
    g = CF( (1+(x-1/3*x**3)**2,(x-1/3*x**3)*(y-1/3*y**3),(x-1/3*x**3)*(z-1/3*z**3),(x-1/3*x**3)*(y-1/3*y**3),1+(y-1/3*y**3)**2,(y-1/3*y**3)*(z-1/3*z**3),(x-1/3*x**3)*(z-1/3*z**3),(y-1/3*y**3)*(z-1/3*z**3),1+(z-1/3*z**3)**2), dims=(3,3) )

    Scalar_ex = -(((2 - 2*z**2)*y**2 + 2*z**2 - 2)*x**12 + ((-24 + 24*z**2)*y**2 - 24*z**2 + 24)*x**10 + ((2 - 2*z**2)*y**6 + (-12 + 12*z**2)*y**4 + (108 - 2*z**6 - 144*z**2 + 12*z**4)*y**2 - 72 - 12*z**4 + 108*z**2 + 2*z**6)*x**8 + ((2 - 2*z**2)*y**8 + (-28 + 28*z**2)*y**6 + (114 - 114*z**2)*y**4 + (468*z**2 + 28*z**6 - 2*z**8 - 198 - 114*z**4)*y**2 - 72 + 114*z**4 + 2*z**8 - 28*z**6 - 198*z**2)*x**6 + ((-12 + 12*z**2)*y**8 + (114 - 114*z**2)*y**6 + (-360 + 360*z**2)*y**4 + (360*z**4 - 114*z**6 + 12*z**8 + 54 - 702*z**2)*y**2 + 594 + 114*z**6 - 12*z**8 - 360*z**4 + 54*z**2)*x**4 + ((2 - 2*z**2)*y**12 + (-24 + 24*z**2)*y**10 + (108 - 2*z**6 - 144*z**2 + 12*z**4)*y**8 + (468*z**2 + 28*z**6 - 2*z**8 - 198 - 114*z**4)*y**6 + (360*z**4 - 114*z**6 + 12*z**8 + 54 - 702*z**2)*y**4 + (486 + 468*z**6 - 2*z**12 + 24*z**10 - 144*z**8 - 702*z**4)*y**2 + 108*z**8 - 324 - 24*z**10 + 2*z**12 - 198*z**6 + 54*z**4 + 486*z**2)*x**2 + (2*z**2 - 2)*y**12 + (-24*z**2 + 24)*y**10 + (-72 - 12*z**4 + 108*z**2 + 2*z**6)*y**8 + (-72 + 114*z**4 + 2*z**8 - 28*z**6 - 198*z**2)*y**6 + (594 + 114*z**6 - 12*z**8 - 360*z**4 + 54*z**2)*y**4 + (108*z**8 - 324 - 24*z**10 + 2*z**12 - 198*z**6 + 54*z**4 + 486*z**2)*y**2 - 486 - 72*z**8 + 24*z**10 - 2*z**12 - 72*z**6 + 594*z**4 - 324*z**2)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)/(x**6 + y**6 + z**6 - 6*x**4 - 6*y**4 - 6*z**4 + 9*x**2 + 9*y**2 + 9*z**2 + 9)**2

    Ricci11 = -(-1 + x**2)*((2*y**2)*z**2 - 2/3*y**2*z**4 + 1/9*y**2*z**6 - 2/3*y**4*z**2 + 1/9*y**6*z**2 + x**2*y**2 - 2/3*x**4*y**2 + 1/9*x**6*y**2 + x**2*z**2 - 2/3*x**4*z**2 + 1/9*x**6*z**2 - 2 - 2/9*x**6 - 1/9*y**6 - 1/9*z**6 + 2/3*y**4 + 4/3*x**4 - 2*x**2 + 2/3*z**4)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6  - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci12 = (1 - 1/3*y**2)*y*(-1 + 1/3*x**2)*x*(x**2 - 1)*(y**2 - 1)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci13 = (1 - 1/3*z**2)*z*(-1 + 1/3*x**2)*x*(x**2 - 1)*(z**2 - 1)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci22 = -(-1 + y**2)*(y**2*z**2 - 2/3*y**4*z**2 + 1/9*y**6*z**2 + x**2*y**2 + (2*x**2)*z**2 - 2/3*x**4*z**2 + 1/9*x**6*z**2 - 2 - 1/9*x**6 - 2/9*y**6 - 1/9*z**6 + 4/3*y**4 - 2*y**2 + 2/3*x**4 - 2/3*x**2*y**4 + 1/9*x**2*y**6 - 2/3*x**2*z**4 + 1/9*x**2*z**6 + 2/3*z**4)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci23 = (1 - 1/3*z**2)*z*(-1 + 1/3*y**2)*y*(y**2 - 1)*(z**2 - 1)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci33 = -(-1 + z**2)*(y**2*z**2 - 2/3*y**2*z**4 + 1/9*y**2*z**6 + (2*x**2)*y**2 - 2/3*x**4*y**2 + 1/9*x**6*y**2 + x**2*z**2 - 2 - 1/9*x**6 - 1/9*y**6 - 2/9*z**6 + 2/3*y**4 + 2/3*x**4 - 2/3*x**2*y**4 + 1/9*x**2*y**6 - 2/3*x**2*z**4 + 1/9*x**2*z**6 + 4/3*z**4 - 2*z**2)/(1 + 1/9*x**6 + 1/9*y**6 + 1/9*z**6 - 2/3*y**4 + y**2 - 2/3*x**4 + x**2 - 2/3*z**4 + z**2)**2
    Ricci_ex = -CF( (Ricci11,Ricci12,Ricci13, Ricci12,Ricci22,Ricci23, Ricci13,Ricci23,Ricci33), dims=(3,3) )

    Einstein_ex = Ricci_ex - 1/2*Scalar_ex*g
    
    fes = HCurlCurl(mesh,order=order)
    gf = GridFunction(fes)
    gf.Set(g)

    assert sqrt(Integrate( InnerProduct(gf.Operator("scalar")-Scalar_ex,gf.Operator("scalar")-Scalar_ex), mesh)) < 1e-7
    assert sqrt(Integrate( InnerProduct(gf.Operator("Ricci")-Ricci_ex,gf.Operator("Ricci")-Ricci_ex), mesh)) < 1e-7
    assert sqrt(Integrate( InnerProduct(gf.Operator("Einstein")-Einstein_ex,gf.Operator("Einstein")-Einstein_ex), mesh)) < 1e-7
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
    Test(mesh=mesh, space=NormalFacetFESpace, order=2, idop = lambda cf : Pn_2d*cf, trace = lambda cf : Pn_2d*cf, vb=VOL, set_dual=[True], facet=True)
    Test(mesh=mesh, space=TangentialFacetFESpace, order=2, idop = lambda cf : Ptau_2d*cf, trace = lambda cf : Ptau_2d*cf, vb=VOL, set_dual=[True], facet=True)
    
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

    mesh = Mesh(unit_square.GenerateMesh(maxh=1,quad_dominated=False))
    TestCurvatureDiffOps2D(mesh,8)
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
    Test(mesh=mesh, space=NormalFacetFESpace, order=2, idop = lambda cf : Pn_3d*cf, trace = lambda cf : Pn_3d*cf, vb=VOL, set_dual=[], facet=True)# dual=True: missing TET element for CalcDual
    
    
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

    mesh = MakeStructured3DMesh(hexes=False, nx=1,ny=1,nz=1, prism=False)
    TestCurvatureDiffOps3D(mesh, order=6)
    
    return


def test_fespaces_surface():
    n_3d = specialcf.normal(3)
    t_3d = specialcf.tangential(3)
    mu_3d = Cross(n_3d,t_3d)
    Ptau_3d = Id(3) - OuterProduct(n_3d,n_3d)
    Pn_3d = OuterProduct(n_3d,n_3d)
    Pt_3d = OuterProduct(t_3d,t_3d)
    Pmu_3d = OuterProduct(mu_3d,mu_3d)

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
    
    #Test(mesh=mesh, space=NormalFacetSurface, order=2, idop = lambda cf : Pmu_3d*cf, diffops=None, vb=BND, set_dual=[True], facet=True) #no dual diffop
    #Test(mesh=mesh, space=HDivSurface, order=0, idop = lambda cf : Ptau_3d*cf, trace = None, diffops=None, vb=BND, set_dual=[False], addorder=0)
    
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
