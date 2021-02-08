git submodule update --init --recursive
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
rm -rf $SRC_DIR
mkdir -p $SRC_DIR
cp -a . $SRC_DIR/
cd $BUILD_DIR

pip3 install pybind11-stubgen==0.5

cmake $SRC_DIR \
      -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_NATIVE_ARCH=OFF \
      -DUSE_CCACHE=ON \
      -DUSE_CGNS=ON \
      -DUSE_UMFPACK=ON \
      -DENABLE_UNIT_TESTS=ON \
      -DCMAKE_OSX_DEPLOYMENT_TARGET=10.14 \
      -DCPACK_PACKAGE_NAME=NGSolve${PACKAGE_NAME_SUFFIX} \
      -DUSE_OCC=ON \
      -DOCC_LIBRARY=/usr/local/opt/opencascade-7.4.0/lib/libTKernel.a \
      -DOCC_INCLUDE_DIR=/usr/local/opt/opencascade-7.4.0/include/opencascade \
      -DOCC_LINK_FREETYPE=ON

make -j5 install
osascript -e 'tell application "Finder" to eject (every disk whose ejectable is true)'
make bundle

