include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# device kernel sources embedded as C++ strings (ngsfmm::code_<kernel>_gpu)
foreach(kernel rotatey shiftz transform p2m l2p m2l_fused)
  add_custom_command(
          OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp
          COMMAND ${CMAKE_COMMAND}
                  -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${kernel}_gpu_kernel.h
                  -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp
                  -DSYMBOL=ngsfmm::code_${kernel}_gpu
                  -P ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd/embed_source.cmake
          DEPENDS ${kernel}_gpu_kernel.h ../ngstd/embed_source.cmake
          COMMENT "Embedding ${kernel}_gpu_kernel.h"
          )
  list(APPEND ngsfmm_gpu_code ${CMAKE_CURRENT_BINARY_DIR}/${kernel}_gpu_code.cpp)
endforeach()

add_library( ngsfmm ${NGS_LIB_TYPE}
        mptools.cpp mptools_device.cpp fmmkernels.cpp fmmoperator.cpp python_fmm.cpp ${ngsfmm_gpu_code}
        )

target_include_directories(ngsfmm PRIVATE ${NETGEN_PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ngstd ${CMAKE_CURRENT_SOURCE_DIR}/../linalg)

target_link_libraries(ngsfmm PUBLIC ngla ngbla ngstd PRIVATE "$<BUILD_INTERFACE:netgen_python>")
install( TARGETS ngsfmm ${ngs_install_dir} )

install( FILES
        mptools.hpp besseljs3d.hpp mptools_device.hpp fmminterface.hpp fmmkernels.hpp fmmoperator.hpp
        normalized_legendre.hpp
        DESTINATION ${NGSOLVE_INSTALL_DIR_INCLUDE}
        COMPONENT ngsolve_devel
       )
