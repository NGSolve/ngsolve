install(FILES
        poisson.py adaptive.py cmagnet.py navierstokes.py elasticity.py
        DESTINATION share/ngsolve/py_tutorials/intro
        COMPONENT ngsolve_devel
       )

install(FILES
        mixed.py hybrid_mixed.py hybrid_dg.py nonlin.py taskmanager.py
        DESTINATION share/ngsolve/py_tutorials/howto
        COMPONENT ngsolve_devel
       )


# install( FILES 
#        ale.pde ale.py all_in_one.py assemble.py bla.py cg.py cube.vol.gz d1_square.pde instat.pde instat.py integrator.py linalg.py mesh.py module1.py square.in2d square.vol
#        DESTINATION share/ngsolve/py_tutorials
#        COMPONENT ngsolve_devel
#       )

