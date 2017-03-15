configure_file(${CMAKE_CURRENT_LIST_DIR}/package_osx_fixup.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/package_osx.cmake IMMEDIATE @ONLY)
add_custom_target(bundle COMMAND ${CMAKE_COMMAND} "-P" "${CMAKE_CURRENT_BINARY_DIR}/package_osx.cmake")
