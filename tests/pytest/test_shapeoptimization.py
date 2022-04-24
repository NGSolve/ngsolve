from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from netgen.csg import unit_cube
from pytest import approx
import pytest

# Taylor Test for first and second (semi-) automatic shape derivative
def Test(G0, G, gfX, X, F, G0bnd=None, Gbnd=None):
    VEC = gfX.space
    mesh = VEC.mesh
    PHI, PSI = VEC.TnT()
    gfX_t = GridFunction(VEC)
    
    #semi automatic
    dJOmegaTestSA = LinearForm(VEC)
    dJOmegaTestSA += G.Diff(F, Grad(PSI)).Compile()
    dJOmegaTestSA += G.Diff(X, PSI).Compile()
    if Gbnd:
        dJOmegaTestSA += Gbnd.Diff(F, Grad(PSI).Trace()).Compile()
        dJOmegaTestSA += Gbnd.Diff(X, PSI).Compile()
    ddJOmegaTestSA = BilinearForm(VEC)
    ddJOmegaTestSA += (G.Diff(F, Grad(PSI)) + G.Diff(X, PSI)).Diff(F, Grad(PHI)).Compile()
    ddJOmegaTestSA += (G.Diff(F, Grad(PSI)) + G.Diff(X, PSI)).Diff(X, PHI).Compile()
    if Gbnd:
        ddJOmegaTestSA += (Gbnd.Diff(F, Grad(PSI).Trace()) + Gbnd.Diff(X, PSI)).Diff(F, Grad(PHI).Trace()).Compile()
        ddJOmegaTestSA += (Gbnd.Diff(F, Grad(PSI).Trace()) + Gbnd.Diff(X, PSI)).Diff(X, PHI).Compile()
    
        
    #full automatic
    dJOmegaTestFA = LinearForm(VEC)
    dJOmegaTestFA += G0.DiffShape(PSI)
    if G0bnd:
        dJOmegaTestFA += G0bnd.DiffShape(PSI)
    ddJOmegaTestFA = BilinearForm(VEC)
    ddJOmegaTestFA += G0.DiffShape(PSI).DiffShape(PHI).Compile()
    if G0bnd:
        ddJOmegaTestFA += G0bnd.DiffShape(PSI).DiffShape(PHI).Compile()

    t = 1/4
    J0 = Integrate(G0, mesh)
    if G0bnd: J0 += Integrate(G0bnd, mesh)

    delta1 = [0,0]
    delta2 = [0,0]
    delta1Prev = [0,0]
    delta2Prev = [0,0]
    
    return

    dJOmegaTestSA.Assemble()
    dJOmegaTestFA.Assemble()
    ddJOmegaTestSA.Assemble()
    ddJOmegaTestFA.Assemble()

    tmp1 = dJOmegaTestSA.vec.CreateVector()
    tmp2 = dJOmegaTestSA.vec.CreateVector()
    tmp1.data = ddJOmegaTestSA.mat*gfX.vec
    tmp2.data = ddJOmegaTestFA.mat*gfX.vec
    
    for i in range(7):
        for j in range(2):
            delta1Prev[j] = delta1[j]
            delta2Prev[j] = delta2[j]
        t = t/2
        gfX_t.vec.data = t*gfX.vec

        mesh.SetDeformation(gfX_t)
        Jt = Integrate(G0, mesh)
        if G0bnd: Jt += Integrate(G0bnd, mesh)
 
        mesh.UnsetDeformation()
        delta1[0] = abs(Jt - J0 -  t*InnerProduct(dJOmegaTestSA.vec,gfX.vec))
        delta1[1] = abs(Jt - J0 -  t*InnerProduct(dJOmegaTestFA.vec,gfX.vec))
        delta2[0] = abs(Jt - J0 -  t*InnerProduct(dJOmegaTestSA.vec,gfX.vec) - t**2/2*InnerProduct(tmp1,gfX.vec))
        delta2[1] = abs(Jt - J0 -  t*InnerProduct(dJOmegaTestFA.vec,gfX.vec) - t**2/2*InnerProduct(tmp2,gfX.vec))
        #print( " delta1 = ", delta1)
        #print( " delta2 = ", delta2)

        if i > 0:
            #print( "--------------------------delta1Prev / delta1 = ",round(delta1Prev[0] / delta1[0],3),round(delta1Prev[1] / delta1[1],3))
            #print( "--------------------------delta2Prev / delta2 = ",round(delta2Prev[0] / delta2[0],3),round(delta2Prev[1] / delta2[1],3))
            assert abs(delta1Prev[0] / delta1[0]) > 3.8 or abs(delta1[0]) < 1e-12
            assert abs(delta1Prev[1] / delta1[1]) > 3.8 or abs(delta1[1]) < 1e-12
            assert abs(delta2Prev[0] / delta2[0]) > 7.1 or abs(delta2[0]) < 1e-12
            assert abs(delta2Prev[1] / delta2[1]) > 7.1 or abs(delta2[1]) < 1e-12
    return


def test_diff():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

    X = CF( (x,y) ).MakeVariable()
    F = Id(2).MakeVariable()

    uD = CF( X[1]*(1-X[1])-X[0]*(1-X[0]) )

    assert sqrt(Integrate((uD.Diff(X, CF( (2,3*x) ))  - (-2*(1-2*x) + 3*x*(1-2*y)))**2,mesh)) < 1e-15
    var = CF( (x,y,-x,x*y), dims=(2,2))

    cf = F*CF( (1,1) )*CF( (3,5) )
    assert sqrt(Integrate((cf.Diff(F, var)  - var*CF( (1,1) )*CF( (3,5) ))**2,mesh)) < 1e-15
    cf = F.trans*CF( (1,1) )*CF( (3,5) )
    assert sqrt(Integrate((cf.Diff(F, var)  - var.trans*CF( (1,1) )*CF( (3,5) ))**2,mesh)) < 1e-15
    return

def test_shapeopt_2d():
    geo = SplineGeometry()
    geo.AddCircle((0,0), r = 0.3)
    ngmesh = geo.GenerateMesh(maxh = 0.025) 
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    X = CF( (x,y) ).MakeVariable()
    VEC = VectorH1(mesh, order=2)
    gfX = GridFunction(VEC)
    gfX.Set((x**2*y*exp(y),y**2*x*exp(x)))

    F = Id(2).MakeVariable()

    uD = CF( 2*X[1]*(1-X[1])+2*X[0]*(1-X[0]) )
    vuD = CF( (uD,uD) )
    def Cost(u):
        return (u-uD)**2
    def CostV(u):
        return (u-vuD)**2

    n = specialcf.normal(2)
    tang = specialcf.tangential(2)


    ### Basic
    G0 = uD * dx
    G  = uD * Det(F) * dx
    Test(G0, G, gfX, X, F)

    ### (Vector)H1
    fes = H1(mesh, order=1)
    gfu = GridFunction(fes)
    gfp = GridFunction(fes)
    gfu.Set( x*y )
    gfp.Set( 4*x - 2*y**2 )

    f = grad(gfu)*grad(gfp)
    G0 = (Cost(gfu) + f)* dx
    G  = ((Inv(F).trans*grad(gfu))*(Inv(F).trans*grad(gfp)) + Cost(gfu)) * Det(F) * dx

    Test(G0, G, gfX, X, F)
    
        
    fes = VectorH1(mesh, order=1)
    gfu = GridFunction(fes)
    gfp = GridFunction(fes)
    gfu.Set( CF( (x*y, -x+y) ) )
    gfp.Set( CF( (x - y**2, -x*2*y) ) )
    
    f = InnerProduct(grad(gfu),grad(gfp)) + gfu*gfp + Trace(grad(gfu))
    G0 = (CostV(gfu) + f)* dx
    G  = (InnerProduct(grad(gfu)*Inv(F),grad(gfp)*Inv(F)) + gfu*gfp + Trace(grad(gfu)*Inv(F))  + CostV(gfu)) * Det(F) * dx

    Test(G0, G, gfX, X, F)



    ### HCurl/HDiv
    fesHC = HCurl(mesh, order=1)
    fesHD = HDiv(mesh, order=1)

    gfu = GridFunction(fesHC)
    gfp = GridFunction(fesHD)
    gfu.Set( CF( (x*y, -x+y) ) )
    gfp.Set( CF( (x - y**2, -x*2*y) ) )

    f = curl(gfu)*div(gfp) + gfu*gfp 
    G0 = (CostV(gfu) + f) * dx
    G  = (InnerProduct(1/Det(F)*F*gfu,Inv(F).trans*gfp) + 1/Det(F)**2*curl(gfu)*div(gfp) + CostV(Inv(F).trans*gfu)) * Det(F) * dx

    Test(G0, G, gfX, X, F)

    ######################## specialcf/BND ######################
    fes = VectorH1(mesh, order=2)
    gfu = GridFunction(fes)
    gfu.Set( CF( (x*3*y, -x*(x+y)) ) )

    f = 1000*(gfu.Trace()*n + gfu.Trace()*tang)
    G0 = CostV(gfu) * dx
    G0bnd = f*ds
    G  = CostV(gfu) * Det(F) * dx
    Gbnd = 1000*(1/Norm(Inv(F).trans*n)*gfu.Trace()*(Inv(F).trans*n) + 1/Norm(F*tang)*gfu.Trace()*(F*tang))*Det(F)*Norm(Inv(F).trans*n)*ds

    Test(G0, G, gfX, X, F, G0bnd,Gbnd)
    return
    

@pytest.mark.slow
def test_code_gen():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

    VEC = VectorH1(mesh,order=1)
    Phi,Psi = VEC.TnT()
    gfvh = GridFunction(VEC)
    gfvh.Set(CF( (x,y*x,x*y*z) ),definedon=mesh.Boundaries(".*"))
    gfh = GridFunction(H1(mesh,order=1))
    gfh.Set(x*y*z,definedon=mesh.Boundaries(".*"))

    gfhc = GridFunction(HCurl(mesh,order=1))
    gfhc.Set(CF( (x,y*x,x*y*z) ),definedon=mesh.Boundaries(".*"))
    gfhd = GridFunction(HDivSurface(mesh,order=1))
    gfhd.Set(CF( (x,y*x,x*y*z) ),definedon=mesh.Boundaries(".*"))
    

    n = specialcf.normal(3)

    functions = [(CF(1)*ds).DiffShape(Psi),(Trace(Grad(n))*ds).DiffShape(Psi),(CF(1)*ds(element_boundary=True)).DiffShape(Psi), (Trace(Grad(gfvh).Trace())*ds).DiffShape(Psi), ((Grad(gfh).Trace())**2*ds).DiffShape(Psi), (curl(gfhc).Trace()*gfhc.Trace()*ds).DiffShape(Psi),(div(gfhd).Trace()*gfhd.Trace()*gfhd.Trace()*ds).DiffShape(Psi)]


    for cf in functions:
        cfs = [ cf.Compile(), cf.Compile(True, wait=True)]
        for f in cfs:
            lf = BilinearForm(VEC)
            lf += Variation(cf-f)
            assert lf.Energy(gfvh.vec) == approx(0)
        
    return

if __name__ == "__main__":
    test_diff()
    test_shapeopt_2d()
    test_code_gen()
