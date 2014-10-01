from os import environ
from sys import path
path.append(environ['NETGENDIR']+'/../lib')

# from libngspy import *

print ("importing ngsolve modules")
from . import ngstd
from . import bla
from . import fem
from . import la
from . import comp
from . import solve
print ("done main init")

del environ
del path
