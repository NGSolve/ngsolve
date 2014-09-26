import sys
try:
    # Linux
    from libmesh import meshing
    from libgeom2d import geom2d
    from libcsg import csg
except:
    # Windows
    from nglib import *

sys.modules['netgen.meshing'] = meshing
sys.modules['netgen.geom2d'] = geom2d
sys.modules['netgen.csg'] = csg
del sys

