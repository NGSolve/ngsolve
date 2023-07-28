# Tests currently only working on Linux
if(NETGEN_USE_PYTHON OR NOT WIN32)

include(CTest)

if(NETGEN_USE_PYTHON )
  if(NETGEN_USE_MPI)
    set(test_driver ngspy ${CMAKE_CURRENT_SOURCE_DIR}/pde_tester.py)
  else()
    set(test_driver ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pde_tester.py)
  endif()
else()
  set(test_driver "${CMAKE_CURRENT_BINARY_DIR}/pde_tester")
  add_executable( pde_tester pde_tester.cpp)
  install( TARGETS pde_tester DESTINATION ${NGSOLVE_INSTALL_DIR_BIN} COMPONENT ngsolve )
  if(NETGEN_USE_GUI)
    target_link_libraries(pde_tester solve  gui netgen_python ${LAPACK_LIBRARIES})
  else()
    target_link_libraries(pde_tester solve netgen_python ${LAPACK_LIBRARIES})
  endif()
endif()

file(COPY ${PROJECT_SOURCE_DIR}/pde_tutorial DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(pde_dir "${CMAKE_CURRENT_BINARY_DIR}/pde_tutorial")
file(COPY ${PROJECT_SOURCE_DIR}/tests/acc_tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(acc_dir "${CMAKE_CURRENT_BINARY_DIR}/acc_tests")



########################################################
# PDE TESTS - SHM
########################################################
set(shm_pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d8_coilshield.pde;d9_hybridDG.pde;d10_DGdoubleglazing.pde;d11_chip_nitsche.pde")

if(USE_UMFPACK)
  if(NETGEN_USE_MPI)
    add_test(NAME umfpack.py COMMAND ngspy umfpack.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    add_test(NAME umfpack.py COMMAND ${NETGEN_PYTHON_EXECUTABLE} umfpack.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  set_tests_properties ( ${py_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
endif(USE_UMFPACK)

foreach(pde_test ${shm_pde_tests})
    # if(EXISTS "${PROJECT_SOURCE_DIR}/tests/pde_tutorial/${pde_test}.appx")
    #     file(READ ${PROJECT_SOURCE_DIR}/tests/pde_tutorial/${pde_test}.appx CONTENTS)
    #     file(APPEND ${pde_dir}/${pde_test} ${CONTENTS})
    # endif()
    add_test(NAME ${pde_test}  COMMAND ${test_driver} "${pde_dir}/${pde_test}" WORKING_DIRECTORY ${pde_dir})
    set_tests_properties ( ${pde_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    set_tests_properties ( ${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( ${pde_test} PROPERTIES LABELS "standard")
endforeach()

add_test(NAME d4_cube_performance.pde COMMAND ${test_driver} "${pde_dir}/d4_cube.pde" "2" WORKING_DIRECTORY ${pde_dir})
set_tests_properties ( d4_cube_performance.pde PROPERTIES LABELS "performance")

set(acc_tests "poisson_circle.pde;poisson_circle_HDG.pde;poisson_circle_HDG_hodc.pde")

foreach(pde_test ${acc_tests})
    add_test(NAME acc_${pde_test}  COMMAND ${test_driver} "${acc_dir}/${pde_test}" WORKING_DIRECTORY ${acc_dir})
    set_tests_properties ( acc_${pde_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    set_tests_properties ( acc_${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( acc_${pde_test} PROPERTIES LABELS "accuracy")
endforeach()

########################################################
# PDE TESTS - MPI
########################################################
# d2_chip.pde  if run, but adaptive refinement does not work
# d8_coilshield.pde
# d10_DGdoubleglazing.pde SEEMS to work, but no idea what it does...
if(NETGEN_USE_MPI)
  set(mpi_pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d9_hybridDG.pde;d11_chip_nitsche.pde")

  foreach(filename ${mpi_pde_tests})
    file(READ ${pde_dir}/${filename} pdefile)
    file(WRITE ${pde_dir}/mpi_${filename} ${pdefile}) 
  endforeach()

  set(mg_to_bddc_files "d2_chip.pde;d4_cube.pde")
  foreach(filename ${mg_to_bddc_files})
    file(READ ${pde_dir}/${filename} pdefile)
    string(REPLACE "multigrid" "bddc" pdefile ${pdefile})
    file(WRITE ${pde_dir}/mpi_${filename} ${pdefile})
  endforeach()

  set(mg_to_direct_files "d5_beam.pde;d6_shaft.pde;d11_chip_nitsche.pde")
  foreach(filename ${mg_to_direct_files})
    file(READ ${pde_dir}/${filename} pdefile)
    string(REPLACE "multigrid" "direct" pdefile ${pdefile})
    file(WRITE ${pde_dir}/mpi_${filename} ${pdefile})
  endforeach()

  set(PROCESSOR_COUNT 4)
  set(test_driver "mpirun" -np ${PROCESSOR_COUNT} --allow-run-as-root ${test_driver})

  foreach(pde_test ${mpi_pde_tests})
    add_test(NAME mpi_${pde_test}  COMMAND ${test_driver} "${pde_dir}/mpi_${pde_test}" WORKING_DIRECTORY ${pde_dir})
    set_tests_properties ( mpi_${pde_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    set_tests_properties ( mpi_${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( mpi_${pde_test} PROPERTIES LABELS "standard")
  endforeach()

  foreach(pde_test ${acc_tests})
    add_test(NAME mpi_acc_${pde_test}  COMMAND ${test_driver} "${acc_dir}/${pde_test}" WORKING_DIRECTORY ${acc_dir})
    set_tests_properties ( mpi_acc_${pde_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    set_tests_properties ( mpi_acc_${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( mpi_acc_${pde_test} PROPERTIES LABELS "accuracy")
  endforeach()

endif()


if(NOT WIN32 AND NOT APPLE)
    add_test(NAME cpu_usage COMMAND top "-n1" "-b" )
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER )
if(BUILD_TYPE_UPPER STREQUAL "DEBUG")
    set_tests_properties ( "d6_shaft.pde" PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    if(NOT NETGEN_USE_MPI)
        set_tests_properties ( "d8_coilshield.pde" PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
    endif(NOT NETGEN_USE_MPI)
endif(BUILD_TYPE_UPPER STREQUAL "DEBUG")

add_test(NAME test_ngscxx  COMMAND ${CMAKE_INSTALL_PREFIX}/${NGSOLVE_INSTALL_DIR_BIN}/${NGSCXX_COMMAND} -c "${CMAKE_CURRENT_SOURCE_DIR}/pde_tester.cpp")

endif(NETGEN_USE_PYTHON OR NOT WIN32)

add_subdirectory(pytest)
add_subdirectory(catch)
add_subdirectory(timings)
