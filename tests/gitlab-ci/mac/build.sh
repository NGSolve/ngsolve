git submodule update --init --recursive
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
rm -rf $SRC_DIR
mkdir -p $SRC_DIR
cp -a . $SRC_DIR/
cd $BUILD_DIR

cmake $SRC_DIR \
      -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_NATIVE_ARCH=OFF \
      -DUSE_CCACHE=ON \
      -DUSE_UMFPACK=ON \
      -DENABLE_UNIT_TESTS=ON \
      -DCMAKE_OSX_DEPLOYMENT_TARGET=10.9 \
      -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
      -DCPACK_PACKAGE_NAME=NGSolve${PACKAGE_NAME_SUFFIX} \

make -j5 install
osascript -e 'tell application "Finder" to eject (every disk whose ejectable is true)'
make bundle

