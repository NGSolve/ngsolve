cd $BUILD_DIR/ngsolve

pip3 install scipy

ctest . --output-on-failure
