install( FILES
        beam.geo chip.vol coilshield.vol cube.vol d3_helmholtz.pde d6_shaft.pde piezo2d40round4.vol.gz square.in2d beam.vol coil.geo coil.vol d1_square.pde d4_cube.pde d7_coil.pde shaft.geo square.vol chip.in2d coilshield.geo cube.geo d2_chip.pde d5_beam.pde d8_coilshield.pde shaft.vol doubleglazing.in2d doubleglazing.vol d10_DGdoubleglazing.pde d9_hybridDG.pde d11_chip_nitsche.pde
        DESTINATION share/ngsolve
       )

include(CTest)

add_test("pde_tutorial_d1" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d1_square.pde")
add_test( "pde_tutorial_d2" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d2_chip.pde")
add_test( "pde_tutorial_d3" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d3_helmholtz.pde")
add_test( "pde_tutorial_d4" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d4_cube.pde")
add_test( "pde_tutorial_d5" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d5_beam.pde")
add_test( "pde_tutorial_d6" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d6_shaft.pde")
add_test( "pde_tutorial_d7" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d7_coil.pde")
add_test( "pde_tutorial_d8" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d8_coilshield.pde")
add_test( "pde_tutorial_d9" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d9_hybridDG.pde")
add_test("pde_tutorial_d10" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d10_DGdoubleglazing.pde")
add_test("pde_tutorial_d11" "${CMAKE_CURRENT_BINARY_DIR}/../tests/ngs_solved" 
  "../solve/ngs" "${PROJECT_SOURCE_DIR}/pde_tutorial/d11_chip_nitsche.pde")
