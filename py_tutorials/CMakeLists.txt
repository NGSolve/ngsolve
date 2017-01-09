install(FILES
        poisson.py adaptive.py cmagnet.py navierstokes.py elasticity.py
        DESTINATION share/ngsolve/py_tutorials/intro
        COMPONENT ngsolve_devel
       )

install(FILES
        mixed.py hybrid_mixed.py hybrid_dg.py nonlin.py taskmanager.py
        pickling.py DG/timeDG-skeleton.py DG/timeDGwave.py
        DG/timeDG.py DG/timeDGlap.py
        DESTINATION share/ngsolve/py_tutorials/howto
        COMPONENT ngsolve_devel
       )

install(FILES
        TensorProduct/TensorProductTools.py TensorProduct/make_seg_mesh.py
        TensorProduct/Bloch_oszillations.py
        TensorProduct/tp_dg_2d.py TensorProduct/tp_dg_2d.py
        DESTINATION share/ngsolve/py_tutorials/TensorProduct
        COMPONENT ngsolve_devel
       )

# install( FILES 
#        ale.pde ale.py all_in_one.py assemble.py bla.py cg.py cube.vol.gz d1_square.pde instat.pde instat.py integrator.py linalg.py mesh.py module1.py square.in2d square.vol
#        DESTINATION share/ngsolve/py_tutorials
#        COMPONENT ngsolve_devel
#       )

include(CTest)
foreach(py_test assemble.py poisson.py adaptive.py cmagnet.py mixed.py hybrid_mixed.py hybrid_dg.py taskmanager.py pickling.py)
  add_test(NAME ${py_test} COMMAND ${NETGEN_PYTHON_EXECUTABLE} ${py_test} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties ( ${py_test} PROPERTIES TIMEOUT 60 )
endforeach()
#  skip those files because the take too long (navierstokes) or need user input (elasticity, nonlin)
#  navierstokes.py elasticity.py nonlin.py 
