from . import comp

from ngsolve import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libsolve.solve import *

if __platform.startswith('win'):
    # Windows
    from ngslib.solve import *

__all__ = [ ]

