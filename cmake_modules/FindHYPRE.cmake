message("hint is ${CMAKE_PREFIX_PATH}")
message("hypre-specifichint is ${HYPRE_HINTS}")

find_path (HYPRE_DIR include/HYPRE.h HINTS "${HYPRE_HINTS}")
if( EXISTS ${HYPRE_DIR}/include/HYPRE.h )
  message( "HYPRE_DIR:  ${HYPRE_DIR}")
  set(HYPRE_FOUND YES)
  # --- includes ---
  set(HYPRE_INCLUDES ${HYPRE_DIR})
  find_path(HYPRE_INCLUDE_DIR HYPRE.h HYPRE_parcsr_ls.h HINTS "${HYPRE_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  list(APPEND HYPRE_INCLUDES ${HYPRE_INCLUDE_DIR})
  message( "HYPRE_INCLUDES = ${HYPRE_INCLUDES}")
  # --- libraries ---
  find_library(LIB_HYPRE_STATIC HYPRE PATHS ${HYPRE_DIR}/lib)
  set(HYPRE_LIBRARIES "${LIB_HYPRE_STATIC}")
  message( "HYPRE_LIBRARIES = ${HYPRE_LIBRARIES}")
else( EXISTS ${HYPRE_DIR}/include/HYPRE.h )
  message ("does not ex: ${HYPRE_DIR}/include/HYPRE.h")
  message( "Could NOT find  HYPRE_DIR:  ${HYPRE_DIR}")
  set(HYPRE_FOUND NO)
endif( EXISTS ${HYPRE_DIR}/include/HYPRE.h )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG HYPRE_LIBRARIES HYPRE_INCLUDES)
  


