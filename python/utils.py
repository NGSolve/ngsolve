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

def L2(mesh, **args):
    """ Create L2 finite element space. """
    return FESpace("l2ho", mesh, **args)

def HCurl(mesh, **args):
    """ Create H(curl) finite element space. """
    return FESpace("hcurlho", mesh, **args)

def HDiv(mesh, **args):
    """ Create H(div) finite element space. """
    return FESpace("hdivho", mesh, **args)



def grad(func):
    if func.derivname != "grad":
        raise Exception("cannot form grad")
    return func.Deriv()

def curl(func):
    if func.derivname != "curl":
        raise Exception("cannot form curl")
    return func.Deriv()

def div(func):
    if func.derivname != "div":
        raise Exception("cannot form div")
    return func.Deriv()




__all__ = ['x', 'y', 'z', 'Laplace', 'Mass', 'Source', 'Neumann', 'H1', 'HCurl', 'HDiv', 'L2', 'grad', 'curl', 'div' ]
