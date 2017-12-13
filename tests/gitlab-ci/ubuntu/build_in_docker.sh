set -e

# Run ssh-agent (inside the build environment)
eval $(ssh-agent -s)
# Add the SSH key stored in SSH_PRIVATE_KEY variable to the agent store
ssh-add <(echo "$SSH_PRIVATE_KEY")
mkdir -p ~/.ssh
[[ -f /.dockerenv ]] && echo -e "Host *\n\tStrictHostKeyChecking no\n\n" > ~/.ssh/config


if [ "$IMAGE_VERSION" == "debug" ]
then
  export CMAKE_CXX_FLAGS="-Og -Wall -Wno-sign-compare -DDebug"
  export CMAKE_BUILD_TYPE="Debug"
else
  export CMAKE_BUILD_TYPE="Release"
fi

if [ "$IMAGE_VERSION" == "avx" ] || [ "$IMAGE_VERSION" == "avx512" ]
then
  export USE_NATIVE_ARCH="ON"
else
  export USE_NATIVE_ARCH="OFF"
fi

cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve \
  -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
  -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
  -DUSE_NATIVE_ARCH=$USE_NATIVE_ARCH \
  -DUSE_OCC=ON \
  -DUSE_CCACHE=ON \
  -DUSE_MKL=ON \
  -DUSE_UMFPACK=ON \
  -DINSTALL_PROFILES=OFF \
  -DMKL_STATIC=ON \
  -DENABLE_UNIT_TESTS=ON \
  -DNG_INSTALL_DIR_LIB=lib/netgen \
  -DNG_INSTALL_DIR_INCLUDE=include/netgen \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DCPACK_PACKAGING_INSTALL_PREFIX=/usr \
  $CMAKE_ARGS

make -j12
make install
make package
cd ngsolve

if [ "$IMAGE_NAME" != "debug" ] && [ "$IMAGE_NAME" != "avx" ]
then
  ## upload built packages to server
  export UPLOAD_DIR=deploy/builds/$CI_PIPELINE_ID/ubuntu/${UBUNTU_VERSION_NAME}_amd64

  apt-get install -y rsync

  if [ "artful" = "$UBUNTU_VERSION_NAME" ]
  then
    echo "build docu"
    make docs
    echo "upload docu"
    rsync -ztrl --del -e ssh \
      --rsync-path="mkdir -p deploy/builds/$CI_PIPELINE_ID/docu/ && rsync" \
      docs/html/* \
      gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$CI_PIPELINE_ID/docu/
  fi
  cd ..

  rsync -ztrl --del -e ssh \
    --rsync-path="mkdir -p $UPLOAD_DIR && rsync" \
    *.deb \
    gitlab-runner@vector.asc.tuwien.ac.at:$UPLOAD_DIR/
fi
