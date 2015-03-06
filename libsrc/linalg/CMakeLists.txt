add_library( la OBJECT
        densemat.cpp polynomial.cpp bfgs.cpp linopt.cpp linsearch.cpp 
        )

set_target_properties(la PROPERTIES POSITION_INDEPENDENT_CODE ON )
