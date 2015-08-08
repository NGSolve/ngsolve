from netgen import __platform
if __platform.startswith('linux') or __platform.startswith('darwin'):
    # Linux or Mac OS X
    from libgeom2d.geom2d import *
#    import libcsgvis.csgvis as csgvis
#    from libcsgvis.csgvis import MouseMove
    from libmesh.meshing import *
if __platform.startswith('win'):
    # Windows
    from nglib.geom2d import *
#    import nglib.csgvis as csgvis
#    from nglib.csgvis import MouseMove
    from nglib.meshing import *




unit_square = SplineGeometry()
pi1 = unit_square.AppendPoint(0,0)
pi2 = unit_square.AppendPoint(1,0)
pi3 = unit_square.AppendPoint(1,1)
pi4 = unit_square.AppendPoint(0,1)
unit_square.Append(["line",pi1,pi2], bc=1)
unit_square.Append(["line",pi2,pi3], bc=2)
unit_square.Append(["line",pi3,pi4], bc=3)
unit_square.Append(["line",pi4,pi1], bc=4)


all = ['SplineGeometry', 'unit_square']




