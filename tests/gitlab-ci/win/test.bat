cd %NETGEN_BUILD_DIR%/ngsolve

pip3 install scipy

ctest -C Release --output-on-failure
