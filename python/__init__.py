"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""
import pkg_resources
import ctypes
import atexit
import os, sys

from . import config

import netgen

if config.is_python_package and sys.platform.startswith('win'):
    os.add_dll_directory(os.path.dirname(netgen.__file__))

if config.is_python_package and config.USE_MKL:
    if sys.platform == 'linux':
        _mkl = pkg_resources.get_distribution('mkl')
        _mkl_rt = _mkl.get_resource_filename('mkl','../../libmkl_rt.so.1')
        ctypes.CDLL(_mkl_rt)
    elif sys.platform.startswith('win'):
        _mkl = pkg_resources.get_distribution('mkl')
        _mkl_rt = os.path.abspath(_mkl.get_resource_filename('mkl','../../Library/bin/mkl_rt.1.dll'))
        ctypes.CDLL(_mkl_rt)

from .ngslib import __version__, ngstd, bla, la, fem, comp, solve

from netgen import Redraw

from pyngcore import BitArray, TaskManager, SetNumThreads, PajeTrace
from .ngstd import Timers, Timer, IntRange
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
    CompressCompound, BoundaryFromVolumeCF, Interpolate, Variation, \
    NumProc, PDE, Integrate, Region, SymbolicLFI, SymbolicBFI, \
    SymbolicEnergy, Mesh, NodeId, ConvertOperator, ORDER_POLICY, VTKOutput, SetHeapSize, \
    SetTestoutFile, ngsglobals, pml, MPI_Init, ContactBoundary, PatchwiseSolve
from .solve import BVP, CalcFlux, Draw, DrawFlux, \
    SetVisualization
from .utils import x, y, z, dx, ds, grad, Grad, curl, div, Deviator, PyId, PyTrace, \
    PyDet, PyCross, PyCof, PyInv, PySym, PySkew, OuterProduct, TimeFunction, Normalize, printmaster
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

# from ngsolve.ngstd import MPIManager
# MPIManager.InitMPI()
mpi_world = MPI_Init()

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

