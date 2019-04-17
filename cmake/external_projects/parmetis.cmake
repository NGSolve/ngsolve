set(PARMETIS_SRC_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/src/project_parmetis)
set(PARMETIS_DIR ${CMAKE_CURRENT_BINARY_DIR}/dependencies/parmetis)

find_file(MPI_H_FILE mpi.h REQUIRED PATHS ${MPI_C_INCLUDE_PATH} ${MPI_C_HEADER_DIR} ${MPI_C_ADDITIONAL_INCLUDE_DIRS} )
get_filename_component(MPI_HDIR ${MPI_H_FILE} DIRECTORY)

ExternalProject_Add(project_parmetis
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dependencies
  URL "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz"
  URL_MD5 f69c479586bf6bb7aff6a9bc0c739628
  DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_dependencies
  PATCH_COMMAND patch -p1 -i ${CMAKE_CURRENT_LIST_DIR}/parmetis.patch
  CMAKE_ARGS
	 -DCMAKE_POSITION_INDEPENDENT_CODE=ON
         -DMPI_INCLUDE_PATH=${MPI_HDIR}
	 -DGKLIB_PATH=${PARMETIS_SRC_DIR}/metis/GKlib
         -DMETIS_PATH=${PARMETIS_SRC_DIR}/metis/
         -DCMAKE_INSTALL_PREFIX=${PARMETIS_DIR}
         -DMETIS_INSTALL=ON
  UPDATE_COMMAND "" # Disable update
  )

set_vars( NGSOLVE_CMAKE_ARGS PARMETIS_DIR )

list(APPEND DEPENDENCIES project_parmetis)
