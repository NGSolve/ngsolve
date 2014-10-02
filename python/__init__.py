from os import environ
from sys import path
from sys import platform as __platform
path.append(environ['NETGENDIR']+'/../lib')

from . import csg
from . import meshing

del environ
del path
