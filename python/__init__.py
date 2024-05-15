"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""
import atexit
import os, sys

from . import config

import netgen

if config.is_python_package and sys.platform.startswith('win'):
    netgen_dir = os.path.dirname(netgen.__file__)
    os.add_dll_directory(netgen_dir)
    os.environ["PATH"] += os.pathsep + netgen_dir

if config.is_python_package and config.USE_MKL:
    import ctypes
    import importlib.metadata
    for f in importlib.metadata.files('intel_openmp'):
        if f.match('*libiomp?.so') or f.match('*libiomp?md.dll'):
            ctypes.CDLL(str(f.locate()))
    for f in importlib.metadata.files('mkl'):
        if f.match('*mkl_rt*'):
            ctypes.CDLL(str(f.locate()))

from .ngslib import __version__, ngstd, bla, la, fem, comp, solve

from netgen import Redraw, TimeFunction

from pyngcore import BitArray, TaskManager, SetNumThreads, PajeTrace, Timers, Timer
from .ngstd import IntRange
from .bla import Matrix, Vector, InnerProduct, Norm
from .la import BaseMatrix, BaseVector, BlockVector, MultiVector, BlockMatrix, \
    CreateVVector, CGSolver, QMRSolver, GMRESSolver, ArnoldiSolver, \
    Projector, DiagonalMatrix, IdentityMatrix, Embedding, PermutationMatrix, \
    ConstEBEMatrix, ParallelMatrix, PARALLEL_STATUS
from .fem import BFI, LFI, CoefficientFunction, Parameter, ParameterC, ET, \
    POINT, SEGM, TRIG, QUAD, TET, PRISM, PYRAMID, HEX, CELL, FACE, EDGE, \
    VERTEX, FACET, ELEMENT, sin, cos, tan, atan, acos, asin, sinh, cosh, \
    exp, log, sqrt, erf, floor, ceil, Conj, atan2, pow, Sym, Skew, Id, Trace, Inv, Det, Cof, Cross, \
    specialcf, BlockBFI, BlockLFI, CompoundBFI, CompoundLFI, BSpline, \
    IntegrationRule, IfPos, VoxelCoefficient, CacheCF, PlaceholderCF
from .comp import VOL, BND, BBND, BBBND, COUPLING_TYPE, ElementId, \
    BilinearForm, LinearForm, GridFunction, Preconditioner, \
    MultiGridPreconditioner, ElementId, FESpace, ProductSpace, H1, HCurl, \
    HDiv, L2, VectorH1, VectorL2, SurfaceL2, TangentialSurfaceL2, HDivDiv, HCurlCurl, HCurlDiv, \
    HDivSurface, HDivDivSurface, FacetFESpace, TangentialFacetFESpace, \
    NormalFacetFESpace, NormalFacetSurface, \
    FacetSurface, VectorSurfaceL2, VectorFacetFESpace, VectorFacetSurface, \
    NodalFESpace, VectorNodalFESpace, H1LumpingFESpace, \
    NumberSpace, Periodic, Discontinuous, Hidden, VectorValued, MatrixValued, Compress, \
    CompressCompound, PlateauFESpace, BoundaryFromVolumeCF, Interpolate, Variation, \
    Integrate, Region, SymbolicLFI, SymbolicBFI, \
    SymbolicEnergy, Mesh, NodeId, ConvertOperator, ORDER_POLICY, VTKOutput, SetHeapSize, \
    SetTestoutFile, ngsglobals, pml, ContactBoundary, PatchwiseSolve, \
    HCurlAMG, APhiHCurlAMG
from .solve import Draw, \
    SetVisualization
from .utils import x, y, z, dx, ds, grad, Grad, curl, div, Deviator, PyId, PyTrace, \
    PyDet, PyCross, PyCof, PyInv, PySym, PySkew, OuterProduct, PrivateSpace, Normalize, printonce

from . import solvers


try:
    from netgen.occ import unit_square, unit_cube
except:
    pass

CF = CoefficientFunction

from math import pi

from builtins import sum as builtin_sum
def sum(iterable, start=None):
    """NGSolve sum function that uses the first element of an iterable as
start argument if no start argument is provided."""
    if start is not None:
        return builtin_sum(iterable, start)
    generator = iter(iterable)
    try:
        first = next(generator)
    except StopIteration:
        return 0
    return builtin_sum(generator, first)

from .timing import Timing

# add flags docu to docstring
def _add_flags_doc(module):
    all_classes = module.__dict__
    for classname in all_classes:
        instance = all_classes[classname]
        try:
            flags_doc = instance.__flags_doc__()
            if instance.__doc__ == None:
                instance.__doc__ = ""
            if not "Keyword arguments can be" in instance.__doc__:
                instance.__doc__ += "\n Keyword arguments can be:\n"
                for name in flags_doc:
                    instance.__doc__ += name + ": " + flags_doc[name] + "\n"
        except AttributeError:
            pass

_add_flags_doc(comp)

# from . import __expr
# BaseVector.expr = property(__expr.VecExpr)
# BaseVector.data = property(__expr.Expr, __expr.expr_data)
# BaseVector.__add__ = __expr.expr_add
# BaseVector.__sub__ = __expr.expr_sub
# BaseVector.__neg__ = __expr.expr_neg
# BaseVector.__rmul__ = __expr.expr_rmul

# BaseMatrix.expr = property(__expr.MatExpr)
# BaseMatrix.data = property(__expr.Expr, __expr.expr_data)
# BaseMatrix.T = property(__expr.TransExpr)
# BaseMatrix.__mul__ = __expr.expr_mul
# BaseMatrix.__rmul__ = __expr.expr_rmul
# BaseMatrix.__neg__ = __expr.expr_neg

def _jupyter_nbextension_paths():
    return [
        {
            "section": "notebook",
            "src": "nbextension/static",
            "dest": "ngsolve_jupyter_widgets",
            "require": "ngsolve_jupyter_widgets/extension",
        }
    ]

atexit.register(solve.__Cleanup)

