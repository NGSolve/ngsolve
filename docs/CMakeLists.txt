if (NETGEN_USE_PYTHON)
  find_package(Sphinx)
  if(SPHINX_FOUND)
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
      ${CMAKE_COMMAND} -E env NETGEN_DOCUMENTATION_RST_FORMAT=1
      ${SPHINX_EXECUTABLE}
      -a -q -b html
      -c "${BINARY_BUILD_DIR}"
      -d "${SPHINX_CACHE_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}"
      "${SPHINX_OUTPUT_DIR}/html"
      COMMENT "Building HTML documentation with Sphinx"
      )

    # not working yet
    # find_package(rst2pdf)
    # if(RST2PDF_FOUND)
    #   add_custom_target(pdf_docs
    # 	${SPHINX_EXECUTABLE}
    # 	-q -b pdf
    # 	-c "${BINARY_BUILD_DIR}"
    # 	-d "${SPHINX_CACHE_DIR}"
    # 	"${CMAKE_CURRENT_SOURCE_DIR}/docs"
    # 	"${SPHINX_OUTPUT_DIR}/pdf"
    # 	COMMENT "Building HTML documentation with Sphinx")
    #   message("rst2pdf found, you can build a pdf documentation with make docs_pdf")
    # else (RST2PDF_FOUND)
    #   message ("rst2pdf not found, if you want to build pdf documentation you have to install it")
    # endif(RSTPDF_FOUND)


    message("Sphinx found, you can build the documentation with make docu")
  else(SPHINX_FOUND)
    message("Sphinx not found, if you want to build the documentation you have to install Sphinx")
  endif(SPHINX_FOUND)

endif(NETGEN_USE_PYTHON)
