configure_file(${CMAKE_CURRENT_SOURCE_DIR}/timings.py ${CMAKE_CURRENT_BINARY_DIR}/timings.py)
add_custom_target(timings ${NETGEN_PYTHON_EXECUTABLE} timings.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
