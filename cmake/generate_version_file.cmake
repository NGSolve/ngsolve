find_package(Git REQUIRED)
if(GIT_FOUND AND EXISTS ${SDIR}/.git)
  execute_process(COMMAND git describe --tags --match "v[0-9]*" --long --dirty WORKING_DIRECTORY ${SDIR} OUTPUT_VARIABLE git_version_string)
else()
  get_filename_component(git_version_string ${SDIR} NAME)
  string(REGEX REPLACE "^ngsolve_(.*)" "\\1" git_version_string "${git_version_string}")
endif()

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" NGSOLVE_VERSION_MAJOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" NGSOLVE_VERSION_MINOR "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" NGSOLVE_VERSION_PATCH "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-([0-9]+).*" "\\1" NGSOLVE_VERSION_TWEAK "${git_version_string}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+\\-[0-9]+\\-([0-9a-z]+).*" "\\1" NGSOLVE_VERSION_HASH "${git_version_string}")

set(NGSOLVE_VERSION_SHORT ${NGSOLVE_VERSION_MAJOR}.${NGSOLVE_VERSION_MINOR}.${NGSOLVE_VERSION_PATCH})
set(NGSOLVE_VERSION_LONG ${NGSOLVE_VERSION_SHORT}-${NGSOLVE_VERSION_TWEAK}-${NGSOLVE_VERSION_HASH})

if(NGSOLVE_VERSION_TWEAK)
  # no release version - nightly build
  set(NGSOLVE_VERSION ${NGSOLVE_VERSION_LONG})
else()
  # TWEAK is 0 -> current version has a tag assigned
  set(NGSOLVE_VERSION ${NGSOLVE_VERSION_SHORT})
endif()

set(NGSOLVE_VERSION_LONG ${NGSOLVE_VERSION_SHORT}-${NGSOLVE_VERSION_TWEAK}-${NGSOLVE_VERSION_HASH})

set(version_file ${BDIR}/ngsolve_version.hpp)
set(new_version_file_string "#define NGSOLVE_VERSION \"${NGSOLVE_VERSION}\"\n")
if(EXISTS ${version_file})
  file(READ ${version_file} old_version_file_string )
  message("old file found")
  if(${old_version_file_string} STREQUAL ${new_version_file_string})
    message("files match, no update")
  else()
    message("files dont match, update")
    file(WRITE ${BDIR}/ngsolve_version.hpp ${new_version_file_string})
  endif()
endif()

