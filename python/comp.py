"""
Finite Element Spaces and Forms
===============================

scalar and vectorial finite element spaces
gridfunctions, bilinear- and linar-forms, preconditioners
"""

from . import fem
from . import la

##from ngsolve import __platform
##if __platform.startswith('linux') or __platform.startswith('darwin'):
##    # Linux or Mac OS X
##    from libngcomp.comp import *

##if __platform.startswith('win'):
##    # Windows
##    from ngslib.comp import *

from ngslib.ngs_docu_flag import SetDocu

if set_docu:
    SetDocu()

from ngslib.comp import *

__all__ = ['BBND', 'BND', 'BilinearForm', 'COUPLING_TYPE', 'CompoundFESpace', 'ElementId', 'BndElementId', 'FESpace', 'GridFunction', 'LinearForm', 'Mesh', 'pml', 'Preconditioner', 'VOL', 'NumProc', 'PDE', 'PyNumProc', 'Integrate', 'IntegrateLF', 'SymbolicLFI', 'SymbolicBFI', 'SymbolicEnergy', 'IntegrateBF', 'VTKOutput', 'ngsglobals']
