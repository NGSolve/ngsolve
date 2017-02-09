# Tests currently only working on Linux
if(NETGEN_USE_PYTHON OR NOT WIN32)

include(CTest)

if(NETGEN_USE_PYTHON)
  set(test_driver ${NETGEN_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pde_tester.py)
else()
  set(test_driver "${CMAKE_CURRENT_BINARY_DIR}/pde_tester")
  add_executable( pde_tester pde_tester.cpp)
  install( TARGETS pde_tester DESTINATION ${ngs_install_dir_bin} COMPONENT ngsolve )
  target_link_libraries(pde_tester solve interface ${TCL_LIBRARY} ${PYTHON_LIBS} ${LAPACK_LIBRARIES})
endif()

file(COPY ${PROJECT_SOURCE_DIR}/pde_tutorial DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(pde_dir "${CMAKE_CURRENT_BINARY_DIR}/pde_tutorial")
file(COPY ${PROJECT_SOURCE_DIR}/tests/acc_tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(acc_dir "${CMAKE_CURRENT_BINARY_DIR}/acc_tests")

set(pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d8_coilshield.pde;d9_hybridDG.pde;d10_DGdoubleglazing.pde;d11_chip_nitsche.pde")

if(USE_UMFPACK)
    add_test(NAME umfpack.py COMMAND ${NETGEN_PYTHON_EXECUTABLE} umfpack.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties ( ${py_test} PROPERTIES TIMEOUT 10 )
endif(USE_UMFPACK)


if(NETGEN_USE_MPI)
# currently d8 and d10 is not working in parallel
    set(pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d9_hybridDG.pde;d11_chip_nitsche.pde")
    if(NOT DEFINED PROCESSOR_COUNT)
        # Unknown:
        set(PROCESSOR_COUNT 0)

        # Linux:
        set(cpuinfo_file "/proc/cpuinfo")
        if(EXISTS "${cpuinfo_file}")
            file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
            list(LENGTH procs PROCESSOR_COUNT)
        endif()

        # Mac:
        if(APPLE)
            find_program(cmd_sys_pro "system_profiler")
            if(cmd_sys_pro)
                execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
                string(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1" PROCESSOR_COUNT "${info}")
            endif()
        endif()

        # Windows:
        if(WIN32)
            set(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
        endif()
    endif()

    set(test_driver "mpirun" -np ${PROCESSOR_COUNT} ${test_driver})

    foreach(pde_test ${pde_tests})
        message("${pde_test}")
        file(READ ${pde_dir}/${pde_test} pdefile)
#no multigrid in parallel implementation
        string(REPLACE "multigrid" "bddc" pdefile ${pdefile})
#no direct solver (MUMPS??)
        string(REPLACE "-solver=direct" "" pdefile ${pdefile})
#segfault at exit if geometry loaded
        string(REPLACE "geometry =" "# geometry =" pdefile ${pdefile})
        file(WRITE ${pde_dir}/${pde_test} ${pdefile})
    endforeach()
    
endif(NETGEN_USE_MPI)

foreach(pde_test ${pde_tests})
    if(EXISTS "${PROJECT_SOURCE_DIR}/tests/pde_tutorial/${pde_test}.appx")
        file(READ ${PROJECT_SOURCE_DIR}/tests/pde_tutorial/${pde_test}.appx CONTENTS)
        file(APPEND ${pde_dir}/${pde_test} ${CONTENTS})
    endif()
    add_test(NAME ${pde_test}  COMMAND ${test_driver} "${pde_dir}/${pde_test}" WORKING_DIRECTORY ${pde_dir})
    set_tests_properties ( ${pde_test} PROPERTIES TIMEOUT 60 )
    set_tests_properties ( ${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( ${pde_test} PROPERTIES LABELS "standard")
endforeach()

    add_test(NAME d4_cube_performance.pde COMMAND ${test_driver} "${pde_dir}/d4_cube.pde" "2" WORKING_DIRECTORY ${pde_dir})
    set_tests_properties ( d4_cube_performance.pde PROPERTIES LABELS "performance")

    set(acc_tests "poisson_circle.pde;poisson_circle_HDG.pde;poisson_circle_HDG_hodc.pde")

foreach(pde_test ${acc_tests})
    add_test(NAME acc_${pde_test}  COMMAND ${test_driver} "${acc_dir}/${pde_test}" WORKING_DIRECTORY ${acc_dir})
    set_tests_properties ( acc_${pde_test} PROPERTIES TIMEOUT 60 )
    set_tests_properties ( acc_${pde_test} PROPERTIES ATTACHED_FILES "${pde_dir}/ng.prof" )
    set_tests_properties ( acc_${pde_test} PROPERTIES LABELS "accuracy")
endforeach()

if(NOT WIN32 AND NOT APPLE)
    add_test(NAME cpu_usage COMMAND top "-n1" "-b" )
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER )
if(BUILD_TYPE_UPPER STREQUAL "DEBUG")
    set_tests_properties ( "d6_shaft.pde" PROPERTIES TIMEOUT 60 )
    if(NOT NETGEN_USE_MPI)
        set_tests_properties ( "d8_coilshield.pde" PROPERTIES TIMEOUT 60 )
    endif(NOT NETGEN_USE_MPI)
endif(BUILD_TYPE_UPPER STREQUAL "DEBUG")

endif(NETGEN_USE_PYTHON OR NOT WIN32)

add_subdirectory(pytest)
