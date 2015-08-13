"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""

from os import environ
from sys import path
path.append(environ['NETGENDIR']+'/../lib')

from sys import platform as __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    path.append(environ['NETGENDIR']+'/../lib')
if __platform.startswith('win'):
    path.append(environ['NETGENDIR'])
    
del environ
del path

__all__ = ['ngstd','bla','fem','la','comp','solve','utils']
 
# ngsmodules = ['ngstd','bla','fem','la','comp','solve']
# variable in import-name ????
# for mod in ngsmodules:
#    from ngsolve.mod import *

# from ngsolve.ngstd import *
# from ngsolve.bla import *
# from ngsolve.la import *
# from ngsolve.fem import *
# from ngsolve.comp import *
# from ngsolve.solve import *
# from ngsolve.utils import *

