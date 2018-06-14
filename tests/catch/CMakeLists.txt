
if(ENABLE_UNIT_TESTS)
if(WIN32)
remove_definitions(-DNGS_EXPORTS)
endif(WIN32)
add_custom_target(unit_tests)

# Build catch_main test object
include_directories(${CATCH_INCLUDE_DIR})
add_library(catch_main STATIC main.cpp)
set_target_properties(catch_main PROPERTIES CXX_STANDARD 11)
add_dependencies(unit_tests catch_main)
add_dependencies(catch_main project_catch)

# ensure the test targets are built before testing
add_test(NAME unit_tests_built COMMAND ${CMAKE_COMMAND} --build . --target unit_tests --config ${CMAKE_BUILD_TYPE} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../.. )

macro(add_unit_test name sources)
  add_executable(test_${name} ${sources} )
  if (WIN32)
    target_link_libraries(test_${name} ngsolve catch_main)
  else(WIN32)
    target_link_libraries(test_${name} ngfem ngstd ngcomp catch_main visual)
  endif(WIN32)

  add_dependencies(unit_tests test_${name})
  add_test(NAME unit_${name} COMMAND test_${name})
  set_tests_properties(unit_${name} PROPERTIES DEPENDS unit_tests_built)
endmacro()

add_unit_test(finiteelement finiteelement.cpp)
add_unit_test(coefficientfunction coefficientfunction.cpp)
add_unit_test(ngblas ngblas.cpp)
file(COPY line.vol square.vol cube.vol DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_unit_test(meshaccess meshaccess.cpp)
endif(ENABLE_UNIT_TESTS)
