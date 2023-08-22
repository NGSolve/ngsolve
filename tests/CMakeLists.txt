# Tests currently only working on Linux
if(NETGEN_USE_PYTHON OR NOT WIN32)

include(CTest)

if(USE_UMFPACK)
  if(NETGEN_USE_MPI)
    add_test(NAME umfpack.py COMMAND ngspy umfpack.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    add_test(NAME umfpack.py COMMAND ${NETGEN_PYTHON_EXECUTABLE} umfpack.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  set_tests_properties ( ${py_test} PROPERTIES TIMEOUT ${NGS_TEST_TIMEOUT} )
endif(USE_UMFPACK)

endif(NETGEN_USE_PYTHON OR NOT WIN32)
add_subdirectory(pytest)
add_subdirectory(catch)
add_subdirectory(timings)
