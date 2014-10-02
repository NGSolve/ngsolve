from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libmesh.meshing import *
if __platform.startswith('win'):
    # Windows
    from nglib.meshing import *
