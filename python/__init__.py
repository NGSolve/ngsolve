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

__all__ = ['ngstd','bla','fem','la','comp','solve'];
 
# print ("importing ngsolve modules")
# from . import ngstd
# from . import bla
# from . import fem
# from . import la
# from . import comp
# from . import solve
# print ("done main init")

