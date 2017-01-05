if(BUILD_NETGEN)
  add_subdirectory(netgen)
  add_dependencies(netgen check_submodules_start)
endif()
