from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libcsg.csg import *
    import libcsgvis.csgvis as csgvis
    from libcsgvis.csgvis import MouseMove
    import libmesh.meshing
#    from libmesh.meshing import *
if __platform.startswith('win'):
    # Windows
    from nglib.csg import *
    import nglib.csgvis as csgvis
    from nglib.csgvis import MouseMove
    from nglib.meshing import *


CSGeometry.VS = csgvis.VS
del csgvis

def VS (obj):
    return obj.VS()

