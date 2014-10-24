from . import bla

##from ngsolve import __platform
##if __platform.startswith('linux') or __platform.startswith('darwin'):
##    # Linux or Mac OS X
##    from libngfem.fem import *

##if __platform.startswith('win'):
##    # Windows
##    from ngslib.fem import *

from ngslib.fem import *

__all__ = ['BFI', 'CoefficientFunction', 'ConstantCF', 'DomainConstantCF', 'ET', 'ElementTransformation', 'FiniteElement', 'H1FE', 'HEX', 'L2FE', 'LFI', 'POINT', 'PRISM', 'PYRAMID', 'PythonCF', 'QUAD', 'SEGM', 'TET', 'TRIG', 'VariableCF' ]


