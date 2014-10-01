try:
    # Linux
    from libngfem.ngfem import *
except:
    # Windows
    from nglib.csg import *
    from nglib.meshing import *
