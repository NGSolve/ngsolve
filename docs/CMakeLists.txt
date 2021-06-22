include (ExternalProject)
if (NETGEN_USE_PYTHON)
  find_package(Sphinx)
  if(SPHINX_FOUND)

    # windows cannot handle symlinks and therefore not build c++ docu
    if(NOT WIN32)
      ExternalProject_Add(
        mylittlengs
        GIT_REPOSITORY https://github.com/NGSolve/mylittlengsolve.git
        TIMEOUT 10
        UPDATE_COMMAND ${GIT_EXECUTABLE} pull
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )
      set_target_properties(mylittlengs PROPERTIES EXCLUDE_FROM_ALL TRUE)
      ExternalProject_Get_Property(mylittlengs source_dir)

      add_custom_target(link_mylittlengs
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${source_dir}" "${CMAKE_CURRENT_SOURCE_DIR}/mylittlengs"
        )
      add_dependencies(link_mylittlengs mylittlengs)
    endif(NOT WIN32)

    # configured documentation tools and intermediate build results
    set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

    # HTML output directory
    set(SPHINX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/conf.py"
      "${BINARY_BUILD_DIR}/conf.py"
      @ONLY)

    add_custom_target(docs_impl
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_RST_FORMAT=1
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_OUT_DIR="${SPHINX_OUTPUT_DIR}/html"
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
      ${SPHINX_EXECUTABLE}
      -a -q -b html
      -c "${BINARY_BUILD_DIR}"
      -d "${SPHINX_CACHE_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}"
      "${SPHINX_OUTPUT_DIR}/html"
      COMMENT "Building HTML documentation with Sphinx"
      )

    add_custom_target(docs_latex_impl
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_RST_FORMAT=1
    	${SPHINX_EXECUTABLE}
    	-q -b latex
    	-c "${BINARY_BUILD_DIR}"
    	-d "${SPHINX_CACHE_DIR}"
    	"${CMAKE_CURRENT_SOURCE_DIR}"
    	"${SPHINX_OUTPUT_DIR}/latex"
    	COMMENT "Building Latex documentation with Sphinx")

      if(USE_MPI)
        add_custom_target(start_ipyparallel
          COMMAND nohup ipcluster start --engines=MPI -n 4 --profile=mpi --debug &
          USES_TERMINAL)
        add_dependencies(docs_impl start_ipyparallel)
        add_dependencies(docs_latex_impl start_ipyparallel)
        add_custom_target(docs
          COMMAND pkill ipcluster
          DEPENDS docs_impl)
        add_custom_target(docs_latex
          COMMAND pkill ipcluster
          DEPENDS docs_latex_impl)
      else(USE_MPI)
        add_custom_target(docs
          DEPENDS docs_impl)
        add_custom_target(docs_latex
          DEPENDS docs_latex_impl)
      endif(USE_MPI)

      if(NOT WIN32)
        add_dependencies(docs_impl link_mylittlengs)
        add_dependencies(docs_latex_impl link_mylittlengs)
      endif(NOT WIN32)

    message("Sphinx found, you can build the documentation with make docu")
  else(SPHINX_FOUND)
    message("Sphinx not found, if you want to build the documentation you have to install Sphinx")
  endif(SPHINX_FOUND)

endif(NETGEN_USE_PYTHON)
