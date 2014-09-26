try:
    # Linux
    from libmesh import meshing
    from libgeom2d import geom2d
    from libcsg import csg
except:
    # Windows
    from nglib import *

