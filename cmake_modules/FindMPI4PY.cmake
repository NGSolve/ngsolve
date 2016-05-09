# - Find MPI4PY
# Find the MPI4PY include directory
#
# This module defines the following variables:
#
#   MPI4PY_FOUND            : True if MPI4PY_INCLUDE_DIR are found
#   MPI4PY_INCLUDE_DIR      : where to find mpi4py/mpi4py.h ..

if(NOT MPI4PY_INCLUDE_DIR)
    execute_process(COMMAND
      "python3" "-c" "import mpi4py; print(mpi4py.get_include())"
      OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
      RESULT_VARIABLE MPI4PY_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(MPI4PY_COMMAND_RESULT)
        message("dependency mpi4py could not be found")
        set(MPI4PY_FOUND FALSE)
    endif(MPI4PY_COMMAND_RESULT)
else(NOT MPI4PY_INCLUDE_DIR)
    set(MPI4PY_FOUND TRUE)
endif(NOT MPI4PY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI4PY DEFAULT_MSG MPI4PY_INCLUDE_DIR)
