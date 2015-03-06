add_definitions(-DNGINTERFACE_EXPORTS)
add_library(gprim OBJECT
        adtree.cpp geom2d.cpp geom3d.cpp geomfuncs.cpp 
        geomtest3d.cpp transform3d.cpp spline.cpp splinegeometry.cpp
        )

set_target_properties( gprim PROPERTIES POSITION_INDEPENDENT_CODE ON )

