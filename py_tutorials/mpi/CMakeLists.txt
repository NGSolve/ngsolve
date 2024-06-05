install(FILES
        mpi_poisson.py mpi_cmagnet.py mpi_navierstokes.py
        mpi_timeDG.py ../__init__.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/ngsolve/demos/mpi
        COMPONENT ngsolve_devel
       )

if(NETGEN_USE_MPI)
     include(CTest)
     foreach(py_test mpi_poisson.py mpi_cmagnet.py mpi_timeDG.py mpi_navierstokes.py)
       add_test(NAME ${py_test} COMMAND mpirun -np 5 --allow-run-as-root ngspy ${py_test} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
       set_tests_properties ( ${py_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
     endforeach()
endif()
