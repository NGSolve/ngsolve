git submodule update --init --recursive
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

cmake $CI_PROJECT_DIR \
      -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_NATIVE_ARCH=OFF \
      -DUSE_CCACHE=ON \
      -DUSE_UMFPACK=ON \
      -DENABLE_UNIT_TESTS=ON \
      -DCPACK_PACKAGE_NAME=NGSolve${PACKAGE_NAME_SUFFIX} \

make -j5 install
osascript -e 'tell application "Finder" to eject (every disk whose ejectable is true)'
make bundle

