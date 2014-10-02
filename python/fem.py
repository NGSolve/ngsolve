from ngsolve import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libngfem.ngfem import *

if __platform.startswith('win'):
    # Windows
    from ngslib.fem import *

