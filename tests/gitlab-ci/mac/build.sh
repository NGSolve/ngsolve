git submodule update --init --recursive
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
rm -rf $SRC_DIR
mkdir -p $SRC_DIR
cp -a . $SRC_DIR/
cd $BUILD_DIR

pip3 install --upgrade pybind11-stubgen

cmake $SRC_DIR \
      -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
      -DCMAKE_BUILD_TYPE=Release \
      -DUSE_NATIVE_ARCH=OFF \
      -DUSE_CCACHE=ON \
      -DUSE_CGNS=ON \
      -DUSE_UMFPACK=ON \
      -DENABLE_UNIT_TESTS=ON \
      -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15 \
      -DCPACK_PACKAGE_NAME=NGSolve${PACKAGE_NAME_SUFFIX} \
      -DPython3_EXECUTABLE=`which python3` \
      -DUSE_OCC=ON

make -j5 install

# eject all mounted .dmg volumes
disks=`diskutil list external | sed -n '/[Ss]cheme/s/.*B *//p'`

if [ "$disks" ]
then
echo "$disks" | while read line ; do
    diskutil unmountDisk /dev/$line
    diskutil eject /dev/$line
  done
fi

make bundle

