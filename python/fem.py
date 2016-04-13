"""
Finite Elements
===============

finite element shape functions, and element-matrix/vector integrators
"""

from . import bla

##from ngsolve import __platform
##if __platform.startswith('linux') or __platform.startswith('darwin'):
##    # Linux or Mac OS X
##    from libngfem.fem import *

##if __platform.startswith('win'):
##    # Windows
##    from ngslib.fem import *

from ngslib.fem import *

__all__ = ['BFI', 'CoefficientFunction',  'DomainConstantCF', 'CoordCF', 'ET', 'ElementTransformation', 'ElementTopology', 'FiniteElement', 'ScalarFE', 'H1FE', 'HEX', 'L2FE', 'LFI', 'POINT', 'PRISM', 'PYRAMID', 'PythonCF', 'QUAD', 'SEGM', 'TET', 'TRIG', 'VERTEX', 'EDGE', 'FACE', 'CELL', 'ELEMENT', 'FACET', 'VariableCF', 'SetPMLParameters', 'sin', 'cos', 'tan', 'atan', 'exp', 'log', 'sqrt', 'Conj', 'specialcf', \
           'BlockBFI', 'BlockLFI', 'CompoundBFI', 'CompoundLFI', 'ConstantCF' \
           ]


