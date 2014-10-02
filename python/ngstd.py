from ngsolve import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libngstd.ngstd import *

if __platform.startswith('win'):
    # Windows
    from ngslib.std import *

