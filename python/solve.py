from ngsolve import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libngspy.ngsolve import *

if __platform.startswith('win'):
    # Windows
    from ngslib.solve import *

__all__ = ['NumProc', 'PDE', 'PyNumProc']

