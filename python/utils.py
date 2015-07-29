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


def H1(mesh, order=1, complex=False, dirichlet=[], dim=1, flags={}):
    """ Create H1 finite element space """
    return FESpace("h1ho", mesh, flags=flags, order=order, complex=complex,
                   dirichlet=dirichlet, dim=dim)

# def H1(mesh, **args):
#     """ Create H1 finite element space """
#     return FESpace("h1ho", mesh, **args)

