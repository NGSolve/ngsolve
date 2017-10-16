from ngsolve.fem import *
from ngsolve.comp import *

x = CoordCF(0)
y = CoordCF(1)
z = CoordCF(2)

def Laplace (coef):
    return BFI("laplace", coef=coef)

def Mass (coef):
    return BFI("mass", coef=coef)


def Source (coef):
    return LFI("source", coef=coef)

def Neumann (coef):
    return LFI("neumann", coef=coef)



def H1(mesh, **args):
    """
    Create H1 finite element space.
    documentation of arguments is available in FESpace.
    """
    fes = FESpace("h1ho", mesh, **args)
    return fes


def VectorH1(mesh, **args):
    """ Create H1 finite element space. """
    fes = FESpace("VectorH1", mesh, **args)
    return fes

def L2(mesh, **args):
    """ Create L2 finite element space. """
    return FESpace("l2ho", mesh, **args)

def SurfaceL2(mesh, **args):
    """ Create L2(boundary) finite element space. """
    return FESpace("l2surf", mesh, **args)

def FacetFESpace(mesh, **args):
    """ Create Facet finite element space. """
    return FESpace("facet", mesh, **args)

def HDivDiv(mesh, **args):
    """ Create H(div-div) finite element space. """
    return FESpace("hdivdiv", mesh, **args)

def NumberSpace(mesh, **args):
    """ Create space of real or complex numbers. """
    return FESpace("number", mesh, **args)


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
    if func.derivname != "curl":
        raise Exception("cannot form curl")
    return func.Deriv()

def div(func):
    if func.derivname == "div":
        return func.Deriv()
    add = func.Operator("div")
    if add:
        return add        
    return func.Deriv()


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

def Det(mat):
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



__all__ = ['x', 'y', 'z', 'Laplace', 'Mass', 'Source', 'Neumann', 'H1', 'VectorH1', 'FacetFESpace', 'L2', 'SurfaceL2', 'HDivDiv', 'NumberSpace', 'grad', 'curl', 'div','Mesh', 'ConstantCF', 'DomainConstantCF', 'Id', 'Trace', 'Det', 'Cross']


