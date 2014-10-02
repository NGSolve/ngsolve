from os import environ
from sys import path
from sys import platform as _platform
path.append(environ['NETGENDIR']+'/../lib')

del environ
del path
