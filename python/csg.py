from libngpy.csg import *
import libngpy.csgvis as csgvis
from libngpy.csgvis import MouseMove
from libngpy.meshing import MeshingParameters


CSGeometry.VS = csgvis.VS
SetBackGroundColor = csgvis.SetBackGroundColor
del csgvis

def VS (obj):
    return obj.VS()



def csg_meshing_func (geom, **args):
    return GenerateMesh (geom, MeshingParameters (**args))

CSGeometry.GenerateMesh = csg_meshing_func


unit_cube = CSGeometry()
unit_cube.Add (OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))

