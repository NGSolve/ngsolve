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
import os
from sys import path

from sys import platform as __platform
if __platform.startswith('linux'):
    path.append(os.path.abspath(os.path.dirname(__file__)+'/../../../'))
if __platform.startswith('win'):
    path.append(os.path.abspath(os.path.dirname(__file__)+'/../../../bin'))
if __platform.startswith('darwin'):
    path.append(os.path.abspath(os.path.dirname(__file__)+'/../../../../../MacOS/'))
    
del environ
del path
del os

#old
# __all__ = ['ngstd','bla','fem','la','comp','solve','utils']


#new:
# from ngsolve.ngstd import *
# from ngsolve.bla import *
# from ngsolve.la import *
# from ngsolve.fem import *
# from ngsolve.comp import *
# from ngsolve.solve import *
# from ngsolve.utils import *

from ngslib import *

def __empty_init(x, *args, **kwargs):
    return

comp.BilinearForm.__new__ = comp.CreateBilinearForm
comp.BilinearForm.__init__ = __empty_init

ngstd.__all__ = ['ArrayD', 'ArrayI', 'BitArray', 'Flags', 'HeapReset', 'IntRange', 'LocalHeap', 'Timers', 'RunWithTaskManager', 'TaskManager', 'SetNumThreads']
bla.__all__ = ['Matrix', 'Vector', 'InnerProduct', 'Norm']
la.__all__ = ['BaseMatrix', 'BaseVector', 'CreateVVector', 'InnerProduct', 'CGSolver', 'QMRSolver', 'GMRESSolver', 'ArnoldiSolver', 'Projector']
fem.__all__ =  ['BFI', 'CoefficientFunction',  'DomainConstantCF',  'Parameter', 'CoordCF', 'ET', 'ElementTransformation', 'ElementTopology', 'FiniteElement', 'ScalarFE', 'H1FE', 'HEX', 'L2FE', 'LFI', 'POINT', 'PRISM', 'PYRAMID', 'QUAD', 'SEGM', 'TET', 'TRIG', 'VERTEX', 'EDGE', 'FACE', 'CELL', 'ELEMENT', 'FACET', 'VariableCF', 'SetPMLParameters', 'sin', 'cos', 'tan', 'atan', 'exp', 'log', 'sqrt', 'Conj', 'atan2', 'pow', 'specialcf', \
           'BlockBFI', 'BlockLFI', 'CompoundBFI', 'CompoundLFI', 'ConstantCF', 'BSpline', \
           'IntegrationRule', 'IfPos' \
           ]
# TODO: fem:'PythonCF' comp:'PyNumProc'
comp.__all__ =  ['BBND','BND', 'BilinearForm', 'COUPLING_TYPE', 'CompoundFESpace', 'ElementId', 'BndElementId', 'FESpace', 'GridFunction', 'LinearForm', 'Mesh', 'NodeId', 'ORDER_POLICY', 'Preconditioner', 'VOL', 'NumProc', 'PDE', 'Integrate', 'SymbolicLFI', 'SymbolicBFI', 'SymbolicEnergy', 'VTKOutput', 'SetHeapSize', 'SetTestoutFile', 'ngsglobals','pml','Periodic']           
solve.__all__ =  ['Redraw', 'BVP', 'CalcFlux', 'Draw', 'DrawFlux', 'SetVisualization']

from ngsolve.ngstd import *
from ngsolve.bla import *
from ngsolve.la import *
from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *
from . import timing

from . import __expr
BaseVector.expr = property(__expr.VecExpr)
BaseVector.data = property(__expr.Expr, __expr.expr_data)
BaseVector.__add__ = __expr.expr_add
BaseVector.__sub__ = __expr.expr_sub
BaseVector.__neg__ = __expr.expr_neg
BaseVector.__rmul__ = __expr.expr_rmul

BaseMatrix.expr = property(__expr.MatExpr)
BaseMatrix.data = property(__expr.Expr, __expr.expr_data)
BaseMatrix.T = property(__expr.TransExpr)
BaseMatrix.__mul__ = __expr.expr_mul
BaseMatrix.__rmul__ = __expr.expr_rmul
BaseMatrix.__neg__ = __expr.expr_neg

Timing = timing.Timing

fem.__doc__ = \
"""Finite Elements
===============

finite element shape functions, and element-matrix/vector integrators
"""


__all__ = ngstd.__all__ + bla.__all__ +la.__all__ + fem.__all__ + comp.__all__ + solve.__all__ + utils.__all__ + ["Timing"]




