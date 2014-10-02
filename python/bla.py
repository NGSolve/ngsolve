from ngsolve import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libngbla.ngbla import *

if __platform.startswith('win'):
    # Windows
    from ngslib.bla import *


__all__ = ['Matrix', 'Vector']

