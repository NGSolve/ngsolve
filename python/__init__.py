from os import environ
from sys import path
from sys import platform as __platform
if __platform.startswith('linux'):
    path.append(environ['NETGENDIR']+'/../lib')
if __platform.startswith('win'):
    path.append(environ['NETGENDIR'])
if __platform.startswith('darwin'):
    path.append(environ['NETGENDIR'])

 
# from libngpy import *
import libngpy
# from libngpy import *

# import libngpy 

# from . import csg
# from . import meshing
# from . import geom2d

del environ
del path
