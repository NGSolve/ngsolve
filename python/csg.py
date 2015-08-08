from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libcsg.csg import *
    import libcsgvis.csgvis as csgvis
    from libcsgvis.csgvis import MouseMove
    from libmesh.meshing import *
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



def csg_meshing_func (geom, **args):
    return GenerateMesh (geom, MeshingParameters (**args))

CSGeometry.GenerateMesh = csg_meshing_func


unit_cube = CSGeometry()
unit_cube.Add (OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))

