from os import environ
from sys import path
from sys import platform as __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    path.append(environ['NETGENDIR']+'/../lib')
if __platform.startswith('win'):
    path.append(environ['NETGENDIR'])


# from libngpy import *
import libngpy 

from . import csg
from . import meshing

del environ
del path
