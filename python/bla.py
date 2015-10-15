"""
Basic Linear Algebra
====================

vectors and matrices

Examples:
x = Vector(5)
x[:] = 2
m = Matrix(3,5)
m[:,:] = 0
m[:,2] = 1

y = m * x
"""


from . import ngstd

##from ngsolve import __platform
##if __platform.startswith('linux') or __platform.startswith('darwin'):
##    # Linux or Mac OS X
##    from libngbla.bla import *

##if __platform.startswith('win'):
##    # Windows
##    from ngslib.bla import *

# from ngslib.bla import *
import ngslib
from ngsolve.bla import *

__all__ = ['Matrix', 'Vector', 'InnerProduct']

