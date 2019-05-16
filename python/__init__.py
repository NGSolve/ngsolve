"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""
import netgen

from ngsolve.ngslib import *
from ngsolve.ngslib import __version__

def TmpRedraw(*args, **kwargs):
    solve._Redraw(*args, **kwargs)
    try:
        import netgen
        import tkinter
        cnt = 0
        while(netgen.gui.win.tk.dooneevent(tkinter._tkinter.DONT_WAIT) and cnt < 100):
            cnt += 1
    except:
        pass

solve.Redraw = TmpRedraw
del TmpRedraw



ngstd.__all__ = ['ArrayD', 'ArrayI', 'BitArray', 'Flags', 'HeapReset', 'IntRange', 'LocalHeap', 'Timers', 'RunWithTaskManager', 'TaskManager', 'SetNumThreads', ]
bla.__all__ = ['Matrix', 'Vector', 'InnerProduct', 'Norm']
la.__all__ = ['BaseMatrix', 'BaseVector', 'BlockVector', 'BlockMatrix', 'CreateVVector', 'InnerProduct', 'CGSolver', 'QMRSolver', 'GMRESSolver', 'ArnoldiSolver', 'Projector', 'IdentityMatrix', 'Embedding', 'PermutationMatrix', 'ConstEBEMatrix']
fem.__all__ =  ['BFI', 'CoefficientFunction', 'Parameter', 'CoordCF', 'ET', 'ElementTransformation', 'ElementTopology', 'FiniteElement', 'MixedFE', 'ScalarFE', 'H1FE', 'HEX', 'L2FE', 'LFI', 'POINT', 'PRISM', 'PYRAMID', 'QUAD', 'SEGM', 'TET', 'TRIG', 'VERTEX', 'EDGE', 'FACE', 'CELL', 'ELEMENT', 'FACET', 'SetPMLParameters', 'sin', 'cos', 'tan', 'atan', 'acos', 'asin', 'sinh', 'cosh', 'exp', 'log', 'sqrt', 'floor', 'ceil', 'Conj', 'atan2', 'pow', 'Sym', 'specialcf', \
           'BlockBFI', 'BlockLFI', 'CompoundBFI', 'CompoundLFI', 'BSpline', \
           'IntegrationRule', 'IfPos' \
           ]
# TODO: fem:'PythonCF' comp:'PyNumProc'
comp.__all__ =  ['BBBND', 'BBND','BND', 'BilinearForm', 'COUPLING_TYPE', 'ElementId', 'BndElementId', 'FESpace','HCurl' , 'GridFunction', 'LinearForm', 'Mesh', 'NodeId', 'ORDER_POLICY', 'Preconditioner', 'MultiGridPreconditioner', 'VOL', 'NumProc', 'PDE', 'Integrate', 'Region', 'SymbolicLFI', 'SymbolicBFI', 'SymbolicEnergy', 'VTKOutput', 'SetHeapSize', 'SetTestoutFile', 'ngsglobals','pml','Periodic','H1','VectorH1','L2','VectorL2','SurfaceL2','HDivDiv','HCurlCurl','HCurlDiv','HDivDivSurface','VectorFacet','FacetFESpace','FacetSurface','HDiv','NumberSpace','HDivSurface','HCurl','Compress','CompressCompound','BoundaryFromVolumeCF', 'MPI_Init']
solve.__all__ =  ['Redraw', 'BVP', 'CalcFlux', 'Draw', 'DrawFlux', 'SetVisualization']

from ngsolve.ngstd import *
from ngsolve.bla import *
from ngsolve.la import *
from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.utils import *
from . import timing


# add flags docu to docstring
all_classes = comp.__dict__
for classname in all_classes:
    instance = all_classes[classname]
    try:
        flags_doc = instance.__flags_doc__()
        if instance.__doc__ == None:
            instance.__doc__ = ""
        instance.__doc__ += "\n Keyword arguments can be:\n"
        for name in flags_doc:
            instance.__doc__ += name + ": " + flags_doc[name] + "\n"
    except AttributeError:
        pass

# from ngsolve.ngstd import MPIManager
# MPIManager.InitMPI()
mpi_world = MPI_Init()

if mpi_world.size>1 and mpi_world.rank==0:
    import ngsolve
    print("importing NGSolve-" + ngsolve.__version__)

from . import __expr
BaseVector.expr = property(__expr.VecExpr)
BaseVector.data = property(__expr.Expr, __expr.expr_data)
BaseVector.__add__ = __expr.expr_add
BaseVector.__sub__ = __expr.expr_sub
BaseVector.__neg__ = __expr.expr_neg
BaseVector.__rmul__ = __expr.expr_rmul

BaseMatrix.expr = property(__expr.MatExpr)
BaseMatrix.data = property(__expr.Expr, __expr.expr_data)
# BaseMatrix.T = property(__expr.TransExpr)
BaseMatrix.__mul__ = __expr.expr_mul
# BaseMatrix.__rmul__ = __expr.expr_rmul
# BaseMatrix.__neg__ = __expr.expr_neg

Timing = timing.Timing

fem.__doc__ = \
"""Finite Elements
===============

finite element shape functions, and element-matrix/vector integrators
"""

# Uncomment this to use patched version of pickle (to regain data pickled somewhere between ~ Feb-Dez 2018)

# register our own memory pickler
# import pickle
# import ngsolve
# pickle._Pickler.dispatch[ngsolve.ngstd._MemoryView] = ngsolve.ngstd._PickleMemory
# pickle._Unpickler.dispatch[b"\xf0"[0]] = ngsolve.ngstd._UnpickleMemory
# # use the python pickler and not cPickle one (cause we can't patch it)
# pickle.Pickler, pickle.Unpickler = pickle._Pickler, pickle._Unpickler
# pickle.dump, pickle.dumps, pickle.load, pickle.loads = pickle._dump, pickle._dumps, pickle._load, pickle._loads


__all__ = ngstd.__all__ + bla.__all__ +la.__all__ + fem.__all__ + comp.__all__ + solve.__all__ + utils.__all__ + ["Timing", "solvers", "meshes", "mpi_world"]




