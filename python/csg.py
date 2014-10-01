try:
    # Linux
    from libcsg.csg import *
    from libmesh.meshing import *
except:
    # Windows
    from nglib.csg import *
    from nglib.meshing import *
