if [ -n "${RUN_TIMINGS}" ];
then
  git submodule update --init --recursive
  mkdir build
  mkdir install
  cd build
  cmake ..
    -DUSE_CCACHE=ON
    -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/install
    -DUSE_NATIVE_ARCH=ON
    -DCMAKE_C_COMPILER=$CMAKE_C_COMPILER
    -DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER
    -DCMAKE_CXX_FLAGS="-ffast-math $CMAKE_FLAGS"
    -DUSE_MKL=ON
    -DMKL_ROOT=/opt/intel/mkl
    -DCMAKE_BUILD_TYPE=Release
  make -j install
  export NETGENDIR=$CI_PROJECT_DIR/install/bin
  export PATH=$CI_PROJECT_DIR/install/bin:$PATH
  export PYTHONPATH=$CI_PROJECT_DIR/install/lib/python3.6/site-packages:.
  export LD_LIBRARY_PATH=$CI_PROJECT_DIR/install/lib:.:$LD_LIBRARY_PATH
  cd ngsolve
  make timings;
fi
