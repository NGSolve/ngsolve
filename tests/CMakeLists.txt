file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved "echo CTEST_FULL_OUTPUT\n\$1 \$2 | tee test_out \nret_code=\$? \ncat test_out\necho \"ret_code=\${ret_code}\"\nFOUND=`grep \"Solved\" test_out`  \nif [ \"\${FOUND}\" != \"\" ] \nthen \n cat ng.prof \n exit \${ret_code} \nelse \n exit 1 \nfi\n")
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/../ngs_solved)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ngs_solved DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
