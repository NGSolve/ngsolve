include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# device kernel sources embedded as C++ strings (ngsbem::code_<kernel>_gpu)
foreach(kernel rotatey shiftz transform p2m l2p m2l_fused)
  add_custom_command(
          OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp
          COMMAND ${CMAKE_COMMAND}
                  -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${kernel}_gpu_kernel.h
                  -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp
                  -DSYMBOL=ngsbem::code_${kernel}_gpu
                  -P ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd/embed_source.cmake
          DEPENDS ${kernel}_gpu_kernel.h ../ngstd/embed_source.cmake
          COMMENT "Embedding ${kernel}_gpu_kernel.h"
          )
  list(APPEND ngsbem_gpu_code ${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp)
endforeach()

add_library( ngsbem ${NGS_LIB_TYPE}
        ${ngsbem_object_libs}
        mptools.cpp potentialtools.cpp python_bem.cpp ngbem.cpp potentialcf.cpp
        integralkernel.cpp fmmoperator.cpp intrules_SauterSchwab.cpp
        analytic_integrals.cpp mptools_device.cpp ${ngsbem_gpu_code}
        )

target_include_directories(ngsbem PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries (ngsbem PUBLIC ngcomp ngfem ngla ngbla ngstd ${MPI_CXX_LIBRARIES} PRIVATE "$<BUILD_INTERFACE:netgen_python>" ${HYPRE_LIBRARIES})
target_link_libraries(ngsbem ${LAPACK_CMAKE_LINK_INTERFACE} "$<BUILD_INTERFACE:ngs_lapack>")
install( TARGETS ngsbem ${ngs_install_dir} )

install( FILES
        mptools.hpp besseljs3d.hpp mptools_device.hpp fmminterface.hpp potentialtools.hpp mp_coefficient.hpp
        kernels.hpp integralkernel.hpp bem_diffops.hpp diffopwithfactor.hpp
        intrules_SauterSchwab.hpp analytic_integrals.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )
