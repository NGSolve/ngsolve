from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.comp import DifferentialSymbol

x = CoordCF(0)
y = CoordCF(1)
z = CoordCF(2)

dx = DifferentialSymbol(VOL)
ds = DifferentialSymbol(BND)


def Laplace (coef):
    return BFI("laplace", coef=coef)

def Mass (coef):
    return BFI("mass", coef=coef)


def Source (coef):
    return LFI("source", coef=coef)

def Neumann (coef):
    return LFI("neumann", coef=coef)



#def H1(mesh, **args):
#    """
#    Create H1 finite element space.
#    documentation of arguments is available in FESpace.
#    """
#    fes = FESpace("h1ho", mesh, **args)
#    return fes


#def VectorH1(mesh, **args):
#    """ Create vector-valued H1 finite element space. """
#    fes = FESpace("VectorH1", mesh, **args)
#    return fes

#def L2(mesh, **args):
#    """ Create L2 finite element space. """
#    return FESpace("l2ho", mesh, **args)

#def VectorL2(mesh, **args):
#    """ Create vector-valued L2 finite element space. """
#    fes = FESpace("VectorL2", mesh, **args)
#    return fes

#def SurfaceL2(mesh, **args):
#    """ Create L2(boundary) finite element space. """
#    return FESpace("l2surf", mesh, **args)

#def FacetFESpace(mesh, **args):
#    """ Create Facet finite element space. """
#    return FESpace("facet", mesh, **args)

#def HDivDiv(mesh, **args):
#    """ Create H(div-div) finite element space. """
#    return FESpace("hdivdiv", mesh, **args)

#def NumberSpace(mesh, **args):
#    """ Create space of real or complex numbers. """
#    return FESpace("number", mesh, **args)


def grad(func):
    if func.derivname == "grad":
        return func.Deriv()
    add = func.Operator("grad")
    if add:
        return add        
    #     if func.derivname != "grad":
    raise Exception("cannot form grad")
    # return func.Deriv()

def curl(func):
    if func.derivname == "curl":
        return func.Deriv()
    add = func.Operator("curl")
    if add:
        return add
    raise Exception("cannot form curl")    

def div(func):
    if func.derivname == "div":
        return func.Deriv()
    try:
        return func.Operator("div")
    except:
        pass
    if func.derivname == "grad" and len(func.dims)==2:  # should check for square
        return Trace(grad(func))
    raise Exception("cannot form div")    


def ConstantCF(val):
    print ("Warning: ConstantCF deprecated, just use CoefficientFunction(val)")
    return CoefficientFunction(val)

def DomainConstantCF(values):
    print ("Warning: DomainConstantCF deprecated, just use CoefficientFunction([values])")
    return CoefficientFunction(values)


def Id(dim):
    return CoefficientFunction( tuple( [1 if i==j else 0 for i in range(dim) for j in range(dim)]), dims=(dim,dim) )

def Trace(mat):
    return sum( [mat[i,i] for i in range(mat.dims[0]) ])

def PyDet(mat):
    if mat.dims[0] == 1:
        return mat[0,0]
    elif mat.dims[0] == 2:
        return mat[0,0]*mat[1,1]-mat[0,1]*mat[1,0]
    elif mat.dims[0] == 3:
        return mat[0,0]*(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1]) \
              +mat[1,0]*(mat[2,1]*mat[0,2]-mat[2,2]*mat[0,1]) \
              +mat[2,0]*(mat[0,1]*mat[1,2]-mat[0,2]*mat[1,1])

def Cross(a,b):
    return CoefficientFunction( (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]) )

def Cof(m):
    if m.dims[0] == 1:
        return CoefficientFunction(1, dims=(1,1))
    elif m.dims[0] == 2:
        return CoefficientFunction( (m[1,1], -m[1,0], -m[0,1], m[0,0]), dims=(2,2) )
    elif m.dims[0] == 3:
        return CoefficientFunction( 
                    ( m[1,1]*m[2,2]-m[2,1]*m[1,2],
                     -m[1,0]*m[2,2]+m[2,0]*m[1,2],
                      m[1,0]*m[2,1]-m[2,0]*m[1,1],
    
                     -m[0,1]*m[2,2]+m[2,1]*m[0,2],
                      m[0,0]*m[2,2]-m[2,0]*m[0,2], 
                     -m[0,0]*m[2,1]+m[2,0]*m[0,1],
    
                      m[0,1]*m[1,2]-m[1,1]*m[0,2], 
                     -m[0,0]*m[1,2]+m[1,0]*m[0,2], 
                      m[0,0]*m[1,1]-m[1,0]*m[0,1] ), dims=(3,3) )

def PyInv(m):
    return 1/Det(m)*Cof(m).trans

def PySym(m):
    return 0.5*(m+m.trans)

def Skew(m):
    return 0.5*(m-m.trans)

def OuterProduct(a, b):
    return CoefficientFunction( tuple([a[i]*b[j] for i in range(a.dim) for j in range(b.dim)]), dims=(a.dim,b.dim) )
## 'L2','H1', 'HDivDiv', 'FacetFESpace', 'VectorL2', 'SurfaceL2', 'NumberSpace', 'VectorH1'
__all__ = ['x', 'y', 'z', 'dx', 'ds', 'Laplace', 'Mass', 'Source', 'Neumann', 'grad', 'curl', 'div','Mesh', 'ConstantCF', 'DomainConstantCF', 'Id', 'Trace', 'PyDet', 'Cross', 'Cof', 'PyInv', 'PySym', 'Skew', 'OuterProduct']


