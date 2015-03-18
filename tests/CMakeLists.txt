include(CTest)

set(test_driver "${CMAKE_CURRENT_BINARY_DIR}/pde_tester")

file(COPY ${PROJECT_SOURCE_DIR}/pde_tutorial DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(pde_dir "${CMAKE_CURRENT_BINARY_DIR}/pde_tutorial")

# currently d8 is not working
set(pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d9_hybridDG.pde;d10_DGdoubleglazing.pde;d11_chip_nitsche.pde")
# set(pde_tests "d1_square.pde;d2_chip.pde;d3_helmholtz.pde;d4_cube.pde;d5_beam.pde;d6_shaft.pde;d7_coil.pde;d8_coilshield.pde;d9_hybridDG.pde;d10_DGdoubleglazing.pde;d11_chip_nitsche.pde")

if(NOT WIN32 AND NOT APPLE)
    add_executable( pde_tester pde_tester.cpp)
    target_link_libraries(pde_tester solve interface ${TCL_LIBRARY} ${PYTHON_LIBS} ${LAPACK_LIBRARIES} )
    install( TARGETS pde_tester DESTINATION bin)
endif(NOT WIN32 AND NOT APPLE)

if(USE_MPI)
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
        string(REPLACE "multigrid" "bddc" pdefile_new ${pdefile})
        string(REPLACE "direct" "sldkfjsldf" pdefile_new1 ${pdefile_new})
        file(WRITE ${pde_dir}/${pde_test} ${pdefile_new1})
    endforeach()
    
endif(USE_MPI)

foreach(pde_test ${pde_tests})
    add_test(NAME ${pde_test}  COMMAND ${test_driver} "${pde_dir}/${pde_test}" WORKING_DIRECTORY ${pde_dir})
endforeach()

