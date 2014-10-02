from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libcsg.csg import *
    from libmesh.meshing import *
if __platform.startswith('win'):
    # Windows
    from nglib.csg import *
    from nglib.meshing import *
