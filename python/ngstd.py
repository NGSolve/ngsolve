##from ngsolve import __platform
##if __platform.startswith('linux') or __platform.startswith('darwin'):
##    # Linux or Mac OS X
##    from libngstd.ngstd import *

##if __platform.startswith('win'):
##    # Windows
##    from ngslib.ngstd import *

# from ngslib.ngstd import *
import ngslib
from ngsolve.ngstd import *

__all__ = ['ArrayD', 'ArrayI', 'BitArray', 'Flags', 'HeapReset', 'IntRange', 'LocalHeap', 'Timers', 'RunWithTaskManager']
