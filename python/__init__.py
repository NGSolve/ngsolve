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
from .ngslib import __version__, ngstd, bla, la, fem, comp, solve

from netgen import Redraw

from pyngcore import BitArray, TaskManager, SetNumThreads
from .ngstd import Timers, Timer, IntRange
from .bla import Matrix, Vector, InnerProduct, Norm
from .la import BaseMatrix, BaseVector, BlockVector, BlockMatrix, \
    CreateVVector, CGSolver, QMRSolver, GMRESSolver, ArnoldiSolver, \
    Projector, IdentityMatrix, Embedding, PermutationMatrix, \
    ConstEBEMatrix, ParallelMatrix, PARALLEL_STATUS
from .fem import BFI, LFI, CoefficientFunction, Parameter, ET, \
    POINT, SEGM, TRIG, QUAD, TET, PRISM, PYRAMID, HEX, CELL, FACE, EDGE, \
    VERTEX, FACET, ELEMENT, sin, cos, tan, atan, acos, asin, sinh, cosh, \
    exp, log, sqrt, floor, ceil, Conj, atan2, pow, Sym, Skew, Id, Trace, Inv, Det, Cof, Cross, \
    specialcf, BlockBFI, BlockLFI, CompoundBFI, CompoundLFI, BSpline, \
    IntegrationRule, IfPos, VoxelCoefficient
from .comp import VOL, BND, BBND, BBBND, COUPLING_TYPE, ElementId, \
    BilinearForm, LinearForm, GridFunction, Preconditioner, \
    MultiGridPreconditioner, ElementId, FESpace, H1, HCurl, \
    HDiv, L2, VectorH1, VectorL2, SurfaceL2, HDivDiv, HCurlCurl, HCurlDiv, \
    HDivSurface, HDivDivSurface, FacetFESpace, TangentialFacetFESpace, \
    NormalFacetFESpace, \
    FacetSurface, VectorSurfaceL2, VectorFacetFESpace, VectorFacetSurface, \
    NodalFESpace, VectorNodalFESpace, \
    NumberSpace, Periodic, Discontinuous, Compress, \
    CompressCompound, BoundaryFromVolumeCF, Variation, \
    NumProc, PDE, Integrate, Region, SymbolicLFI, SymbolicBFI, \
    SymbolicEnergy, Mesh, NodeId, ORDER_POLICY, VTKOutput, SetHeapSize, \
    SetTestoutFile, ngsglobals, pml, MPI_Init, ContactBoundary, PatchwiseSolve
from .solve import BVP, CalcFlux, Draw, DrawFlux, \
    SetVisualization
from .utils import x, y, z, dx, ds, grad, Grad, curl, div, PyId, PyTrace, \
    PyDet, PyCross, PyCof, PyInv, PySym, PySkew, OuterProduct, TimeFunction, Normalize
from . import solvers


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


