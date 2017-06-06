"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""

from ngsolve.ngslib import *


def __empty_init(x, *args, **kwargs):
    return

def __hcurl_new(class_t, *args,**kwargs):
    return comp.CreateFESpace(class_t,"hcurlho",*args,**kwargs)

def __mynew(thisclass, creatorfunction):
    pybind_constructor = thisclass.__new__
    def foo(class_t, *args,**kwargs):
        if class_t is thisclass:
            return  creatorfunction(class_t, *args, **kwargs)
        else:
            return pybind_constructor(class_t,*args,**kwargs)
    return foo

# assign creator functions to __new__
comp.BilinearForm.__new__ = comp.CreateBilinearForm
comp.BilinearForm.__init__ = __empty_init

fem.CoefficientFunction.__new__ = __mynew(fem.CoefficientFunction,fem.CreateCoefficientFunction)
fem.CoefficientFunction.__init__ = __empty_init

comp.GridFunction.__new__ = comp.CreateGridFunction
comp.GridFunction.__init__ = __empty_init

comp.LinearForm.__new__ = comp.CreateLinearForm
comp.LinearForm.__init__ = __empty_init

comp.PDE.__new__ = comp.CreatePDE
comp.PDE.__init__ = __empty_init

comp.VTKOutput.__new__ = comp.CreateVTKOutput
comp.VTKOutput.__init__ = __empty_init

fem.ElementTransformation.__new__ = fem.CreateElementTransformation
fem.ElementTransformation.__init__ = __empty_init

fem.BFI.__new__ = fem.CreateBilinearFormIntegrator
fem.BFI.__init__ = __empty_init

fem.LFI.__new__ = fem.CreateLinearFormIntegrator
fem.LFI.__init__ = __empty_init

comp.FESpace.__new__ = comp.CreateFESpace
comp.FESpace.__init__ = __empty_init

comp.HCurl.__new__ = __hcurl_new
comp.HCurl.__init__ = __empty_init

comp.Periodic.__new__ = comp.CreatePeriodicFESpace
comp.Periodic.__init__ = __empty_init

def TmpRedraw(*args, **kwargs):
    solve._Redraw(*args, **kwargs)
    try:
        import netgen
        import tkinter
        while(netgen.gui.win.tk.dooneevent(tkinter._tkinter.DONT_WAIT)):
            pass
    except:
        pass

solve.Redraw = TmpRedraw
del TmpRedraw

ngstd.__all__ = ['ArrayD', 'ArrayI', 'BitArray', 'Flags', 'HeapReset', 'IntRange', 'LocalHeap', 'Timers', 'RunWithTaskManager', 'TaskManager', 'SetNumThreads']
bla.__all__ = ['Matrix', 'Vector', 'InnerProduct', 'Norm']
la.__all__ = ['BaseMatrix', 'BaseVector', 'CreateVVector', 'InnerProduct', 'CGSolver', 'QMRSolver', 'GMRESSolver', 'ArnoldiSolver', 'Projector']
fem.__all__ =  ['BFI', 'CoefficientFunction', 'Parameter', 'CoordCF', 'ET', 'ElementTransformation', 'ElementTopology', 'FiniteElement', 'ScalarFE', 'H1FE', 'HEX', 'L2FE', 'LFI', 'POINT', 'PRISM', 'PYRAMID', 'QUAD', 'SEGM', 'TET', 'TRIG', 'VERTEX', 'EDGE', 'FACE', 'CELL', 'ELEMENT', 'FACET', 'SetPMLParameters', 'sin', 'cos', 'tan', 'atan', 'exp', 'log', 'sqrt', 'Conj', 'atan2', 'pow', 'specialcf', \
           'BlockBFI', 'BlockLFI', 'CompoundBFI', 'CompoundLFI', 'BSpline', \
           'IntegrationRule', 'IfPos' \
           ]
# TODO: fem:'PythonCF' comp:'PyNumProc'
comp.__all__ =  ['BBND','BND', 'BilinearForm', 'COUPLING_TYPE', 'ElementId', 'BndElementId', 'FESpace','HCurl' , 'GridFunction', 'LinearForm', 'Mesh', 'NodeId', 'ORDER_POLICY', 'Preconditioner', 'VOL', 'NumProc', 'PDE', 'Integrate', 'SymbolicLFI', 'SymbolicBFI', 'SymbolicEnergy', 'VTKOutput', 'SetHeapSize', 'SetTestoutFile', 'ngsglobals','pml','Periodic']           
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




