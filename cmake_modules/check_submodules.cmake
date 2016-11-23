execute_process(COMMAND git submodule status --recursive WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" OUTPUT_VARIABLE git_status_output)
string(REPLACE "\n" ";" git_status_output "${git_status_output}")
foreach( a ${git_status_output})
  if(NOT ${a} MATCHES " [a-f,0-9]* external_dependencies[.]*")
    message(WARNING
"*****************************************************************
      WARNING: The git submodules are out of sync! Please run
      git submodule update --init --recursive
      in your source directory
*****************************************************************")
  endif()
endforeach()

