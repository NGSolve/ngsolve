if(NOT BDIR)
  set(BDIR ${CMAKE_CURRENT_BINARY_DIR})
endif()

find_package(Git REQUIRED)

if(GIT_FOUND AND EXISTS ${CMAKE_CURRENT_LIST_DIR}/../.git)
    execute_process(COMMAND git describe --tags --match "v[0-9]*" --long --dirty
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE git_version_string
        RESULT_VARIABLE status
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/../version.txt)
    file(READ ${CMAKE_CURRENT_LIST_DIR}/../version.txt git_version_string )
  else()
    get_filename_component(git_version_string ${CMAKE_CURRENT_LIST_DIR}/.. NAME)
    string(REGEX REPLACE "^ngsolve_(.*)" "\\1" git_version_string "${git_version_string}")
  endif()
endif()
string(STRIP ${git_version_string} git_version_string)

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

if(NOT NGSOLVE_VERSION_PYTHON)
    # Derive a PEP 440 compliant python package version from git, matching
    # tests/utils.py:get_version() (used for the wheel metadata version).
    if(NGSOLVE_VERSION_TWEAK)
        # commits after the last tag -> post release
        set(NGSOLVE_VERSION_PYTHON "${NGSOLVE_VERSION_SHORT}.post${NGSOLVE_VERSION_TWEAK}")
        set(_ngs_dev_build TRUE)
        if(DEFINED ENV{NG_NO_DEV_PIP_VERSION})
            set(_ngs_dev_build FALSE)
        endif()
        if("$ENV{CI_COMMIT_REF_NAME}" STREQUAL "release")
            set(_ngs_dev_build FALSE)
        endif()
        if(_ngs_dev_build)
            set(NGSOLVE_VERSION_PYTHON "${NGSOLVE_VERSION_PYTHON}.dev0")
        endif()
    else()
        # current commit is tagged -> clean release version
        set(NGSOLVE_VERSION_PYTHON "${NGSOLVE_VERSION_SHORT}")
    endif()
endif()

set(version_file ${BDIR}/ngsolve_version.hpp)
set(new_version_file_string "#define NGSOLVE_VERSION \"${NGSOLVE_VERSION}\"\n")
if(EXISTS ${version_file})
  file(READ ${version_file} old_version_file_string )
  if(${old_version_file_string} STREQUAL ${new_version_file_string})
  else()
    file(WRITE ${BDIR}/ngsolve_version.hpp ${new_version_file_string})
  endif()
else()
    file(WRITE ${BDIR}/ngsolve_version.hpp ${new_version_file_string})
endif()

