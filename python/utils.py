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

def HCurl(mesh, **args):
   """ Create H(curl) finite element space. """
   return FESpace("hcurlho", mesh, **args)

def HDiv(mesh, **args):
    """ Create H(div) finite element space. """
    return FESpace("hdivho", mesh, **args)

def FacetFESpace(mesh, **args):
    """ Create Facet finite element space. """
    return FESpace("facet", mesh, **args)

def HDivDiv(mesh, **args):
    """ Create H(div-div) finite element space. """
    return FESpace("hdivdiv", mesh, **args)


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


__all__ = ['x', 'y', 'z', 'Laplace', 'Mass', 'Source', 'Neumann', 'H1', 'VectorH1', 'FacetFESpace', 'HCurl', 'HDiv', 'L2', 'SurfaceL2', 'HDivDiv', 'grad', 'curl', 'div','NgsPickler', 'NgsUnpickler', 'Mesh', 'ConstantCF', 'DomainConstantCF' ]


