if(USE_PYTHON)
    install (FILES
            bla.py comp.py __expr.py fem.py internal.py
            __init__.py la.py ngstd.py solve.py utils.py solvers.py
            DESTINATION ${PYTHON_PACKAGES_INSTALL_DIR}/ngsolve
            COMPONENT ngsolve
            )
endif(USE_PYTHON)
