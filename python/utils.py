from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.comp import HCurlFunctionsWrap

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

class HCurl(HCurlFunctionsWrap):
    def __init__(self,mesh,*args,**kwargs):
        # fix for pickling
        #if mesh=="hcurlho":
        #    FESpace.__init__(self,mesh,*args,**kwargs)
        #else:
        FESpace.__init__(self,"hcurlho",mesh,*args,**kwargs)

#def HCurl(mesh, **args):
#    """ Create H(curl) finite element space. """
#    return FESpace("hcurlho", mesh, **args)

def HDiv(mesh, **args):
    """ Create H(div) finite element space. """
    return FESpace("hdivho", mesh, **args)

def FacetFESpace(mesh, **args):
    """ Create Facet finite element space. """
    return FESpace("facet", mesh, **args)

def HDivDiv(mesh, **args):
    """ Create H(div-div) finite element space. """
    return FESpace("hdivdiv", mesh, **args)


def PeriodicH1(*args,**kwargs):
    """
Periodic H1
===========

Generator function for periodic H1 space. The mesh needs to be periodic, i.e. use the function CSGeometry.PeriodicSurfaces(master,slave) to create a periodic geometry. With this Netgen creates a periodic mesh, the returned function space maps the slave dofs to the master dofs.

    """
    return FESpace("perH1",*args,**kwargs)

def PeriodicHCurl(*args,**kwargs):
    """
Periodic HCurl
===========

Generator function for periodic HCurl space. The mesh needs to be periodic, i.e. use the function CSGeometry.PeriodicSurfaces(master,slave) to create a periodic geometry. With this Netgen creates a periodic mesh, the returned function space maps the slave dofs to the master dofs.

    """
    return FESpace("perHCurl",*args,**kwargs)

def PeriodicHDiv(*args,**kwargs):
    """
Periodic HDiv
===========

Generator function for periodic HDiv space. The mesh needs to be periodic, i.e. use the function CSGeometry.PeriodicSurfaces(master,slave) to create a periodic geometry. With this Netgen creates a periodic mesh, the returned function space maps the slave dofs to the master dofs.

    """
    return FESpace("perHDiv",*args,**kwargs)



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


import pickle

def NgsPickler(*args, **kargs):
    pickler = pickle.Pickler(*args, **kargs)
    dumped_pids = []

    def my_persistent_id(obj):
        try:
            pid = obj.__ngsid__()
            if pid in dumped_pids:
                return dumped_pids.index(pid)
            else:
                dumped_pids.append(pid)
                obj.__persistent_id__ = dumped_pids.index(pid)
                return obj
        except:
            return None

    pickler.persistent_id = my_persistent_id
    return pickler

def NgsUnpickler(*args, **kargs):
    unpickler = pickle.Unpickler(*args, **kargs)
    loaded_pids = {}

    def my_persistent_load(pid):
        if hasattr(pid,'__ngsid__'):
            loaded_pids[pid.__persistent_id__] = pid
            del pid.__persistent_id__
            return pid
        else:
            return loaded_pids[pid]

    unpickler.persistent_load = my_persistent_load
    return unpickler


__all__ = ['x', 'y', 'z', 'Laplace', 'Mass', 'Source', 'Neumann', 'H1', 'FacetFESpace', 'HCurl', 'HDiv', 'L2', 'HDivDiv', 'PeriodicH1', 'PeriodicHCurl', 'PeriodicHDiv', 'grad', 'curl', 'div','NgsPickler', 'NgsUnpickler', 'Mesh' ]
