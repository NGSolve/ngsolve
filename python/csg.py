from netgen import _platform
if _platform.startswith('linux') or _platform.startswith('darwin'):
    # Linux or Mac OS X
    from libcsg.csg import *
    from libmesh.meshing import *
if _platform.startswith('win'):
    # Windows
    from nglib.csg import *
    from nglib.meshing import *
