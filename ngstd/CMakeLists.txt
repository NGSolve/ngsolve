include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(ngstd_sources
        exception.cpp table.cpp bitarray.cpp flags.cpp 
        symboltable.cpp blockalloc.cpp evalfunc.cpp templates.cpp  
        localheap.cpp stringops.cpp profiler.cpp archive.cpp
        cuda_ngstd.cpp python_ngstd.cpp taskmanager.cpp
        paje_interface.cpp bspline.cpp
        )

if(NOT WIN32)
    set(ngstd_sources ${ngstd_sources} sockets.cpp)
endif(NOT WIN32)

add_library( ngstd ${NGS_LIB_TYPE} ${ngstd_sources} )
target_compile_definitions(ngstd PUBLIC ${NGSOLVE_COMPILE_DEFINITIONS})
target_include_directories(ngstd PUBLIC ${NETGEN_INCLUDE_DIRS})

add_dependencies( ngstd generate_version_file )

if(NOT WIN32)
    target_link_libraries(ngstd ${MPI_CXX_LIBRARIES} ${NETGEN_PYTHON_LIBRARIES} ${NUMA_LIB})
    install( TARGETS ngstd ${ngs_install_dir} )
endif(NOT WIN32)

install( FILES
        array.hpp autodiffdiff.hpp autodiff.hpp autoptr.hpp    
        bessel.hpp bitarray.hpp blockalloc.hpp evalfunc.hpp exception.hpp     
        flags.hpp hashtable.hpp localheap.hpp memusage.hpp ngstd.hpp	      
        profiler.hpp stringops.hpp symboltable.hpp table.hpp templates.hpp    
        parthreads.hpp statushandler.hpp ngsstream.hpp mpiwrapper.hpp	      
        polorder.hpp archive.hpp archive_base.hpp sockets.hpp cuda_ngstd.hpp  
        mycomplex.hpp tuple.hpp paje_interface.hpp python_ngstd.hpp ngs_utils.hpp
        taskmanager.hpp bspline.hpp simd.hpp simd_complex.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )
