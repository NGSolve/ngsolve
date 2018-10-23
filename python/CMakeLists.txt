if(NETGEN_USE_PYTHON)
    install (FILES
            __expr.py internal.py __console.py
            __init__.py utils.py solvers.py eigenvalues.py meshes.py
            krylovspace.py nonlinearsolvers.py timing.py TensorProductTools.py
            DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve
            COMPONENT ngsolve
            )
endif(NETGEN_USE_PYTHON)
