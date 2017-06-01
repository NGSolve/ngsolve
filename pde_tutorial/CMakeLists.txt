install( FILES
        beam.geo chip.vol coilshield.vol cube.vol d3_helmholtz.pde d6_shaft.pde piezo2d40round4.vol.gz square.in2d beam.vol coil.geo coil.vol d1_square.pde d4_cube.pde d7_coil.pde shaft.geo square.vol chip.in2d coilshield.geo cube.geo d2_chip.pde d5_beam.pde d8_coilshield.pde shaft.vol doubleglazing.in2d doubleglazing.vol d10_DGdoubleglazing.pde d9_hybridDG.pde d11_chip_nitsche.pde
        DESTINATION ${NGSOLVE_INSTALL_DIR_RES}/ngsolve
        COMPONENT ngsolve_tutorial
       )

