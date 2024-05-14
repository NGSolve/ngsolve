install(FILES
        poisson.py adaptive.py cmagnet.py navierstokes.py
        elasticity.py poisson.ipynb
        DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsolve/py_tutorials/intro
        COMPONENT ngsolve_devel
       )

install(FILES
        mixed.py hybrid_dg.py nonlin.py taskmanager.py
        pickling.py pml.py DG/timeDG-skeleton.py DG/timeDGwave.py
        DG/timeDG.py DG/timeDGlap.py hdivdiv/hhj.py hdivdiv/tdnns.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsolve/py_tutorials/howto
        COMPONENT ngsolve_devel
       )

install(FILES
        TensorProduct/tp_dg_1d_1d.py
        TensorProduct/tp_dg_2d_1d.py TensorProduct/tp_dg_1d_2d.py
        TensorProduct/tp_dg_2d_2d.py
        DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsolve/py_tutorials/TensorProduct
        COMPONENT ngsolve_devel
       )

include(CTest)
foreach(py_test assemble.py bla.py poisson.py adaptive.py cmagnet.py mixed.py hybrid_dg.py taskmanager.py compound.py pickling.py)
  add_test(NAME ${py_test} COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${py_test} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties ( ${py_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
endforeach()

add_subdirectory(mpi)

#  skip those files because the take too long (navierstokes) or need user input (elasticity, nonlin)
#  navierstokes.py elasticity.py nonlin.py 
