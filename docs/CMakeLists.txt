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

    add_custom_target(docs
      COMMAND
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

    add_custom_target(docs_latex
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_RST_FORMAT=1
    	${SPHINX_EXECUTABLE}
    	-q -b latex
    	-c "${BINARY_BUILD_DIR}"
    	-d "${SPHINX_CACHE_DIR}"
    	"${CMAKE_CURRENT_SOURCE_DIR}"
    	"${SPHINX_OUTPUT_DIR}/latex"
    	COMMENT "Building Latex documentation with Sphinx")

    # technotes/<nr>_<name>/<name>.tex -> html/technotes/<nr>_<name>.pdf and .html
    find_program(LATEXMK_EXECUTABLE latexmk)
    find_program(PANDOC_EXECUTABLE pandoc)
    if(NOT PANDOC_EXECUTABLE)
      message("pandoc not found, technotes html pages are not built")
    endif()
    if(LATEXMK_EXECUTABLE)
      add_custom_target(technotes
        COMMAND ${CMAKE_COMMAND} -E make_directory "${SPHINX_OUTPUT_DIR}/html/technotes")
      file(GLOB technote_dirs "${CMAKE_CURRENT_SOURCE_DIR}/technotes/*")
      foreach(dir ${technote_dirs})
        get_filename_component(note "${dir}" NAME)
        string(REGEX REPLACE "^[0-9]+-[0-9]+_" "" name "${note}")
        if(EXISTS "${dir}/${name}.tex")
          set(outdir "${CMAKE_CURRENT_BINARY_DIR}/technotes/${note}")
          add_custom_command(TARGET technotes POST_BUILD
            COMMAND ${LATEXMK_EXECUTABLE} -pdf -quiet -interaction=nonstopmode -outdir=${outdir} ${name}.tex
            COMMAND ${CMAKE_COMMAND} -E copy "${outdir}/${name}.pdf" "${SPHINX_OUTPUT_DIR}/html/technotes/${note}.pdf"
            WORKING_DIRECTORY "${dir}"
            COMMENT "Building technote ${note}")
          if(PANDOC_EXECUTABLE)
            add_custom_command(TARGET technotes POST_BUILD
              COMMAND ${NETGEN_PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/technotes/tex2html.py"
                "${dir}/${name}.tex" "${SPHINX_OUTPUT_DIR}/html/technotes/${note}.html"
              COMMENT "Building technote ${note} (html)")
          endif()
        endif()
      endforeach()
      add_dependencies(docs technotes)
    else()
      message("latexmk not found, technotes PDFs are not built")
    endif()

      if(NOT WIN32)
        add_dependencies(docs link_mylittlengs)
        add_dependencies(docs_latex link_mylittlengs)
      endif(NOT WIN32)

    message("Sphinx found, you can build the documentation with make docu")
  else(SPHINX_FOUND)
    message("Sphinx not found, if you want to build the documentation you have to install Sphinx")
  endif(SPHINX_FOUND)

endif(NETGEN_USE_PYTHON)
