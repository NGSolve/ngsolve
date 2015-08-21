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
    return FESpace("h1ho", mesh, **args)

def L2(mesh, **args):
    """ Create L2 finite element space. """
    return FESpace("l2ho", mesh, **args)

def HCurl(mesh, **args):
    """ Create H(curl) finite element space. """
    return FESpace("hcurlho", mesh, **args)

def HDiv(mesh, **args):
    """ Create H(div) finite element space. """
    return FESpace("hdivho", mesh, **args)

