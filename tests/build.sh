set -e

if [ "$UBUNTU_VERSION_NAME" == "debug" ]
then
  /root/src/ngsolve/tests/build_debug.sh
  exit
fi

# Run ssh-agent (inside the build environment)
eval $(ssh-agent -s)
# Add the SSH key stored in SSH_PRIVATE_KEY variable to the agent store
ssh-add <(echo "$SSH_PRIVATE_KEY")
mkdir -p ~/.ssh
[[ -f /.dockerenv ]] && echo -e "Host *\n\tStrictHostKeyChecking no\n\n" > ~/.ssh/config

export UPLOAD_DIR=deploy/builds/$CI_PIPELINE_ID/ubuntu/${UBUNTU_VERSION_NAME}_amd64
rsync -ztrl --del -e ssh \
  --rsync-path="mkdir -p $UPLOAD_DIR && rsync" \
  *.deb \
  gitlab-runner@vector.asc.tuwien.ac.at:$UPLOAD_DIR/

cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve \
  -DUSE_NATIVE_ARCH=OFF \
  -DUSE_OCC=ON \
  -DUSE_CCACHE=ON \
  -DUSE_MKL=ON \
  -DUSE_UMFPACK=ON \
  -DINSTALL_PROFILES=OFF \
  -DMKL_STATIC=ON \
  -DENABLE_UNIT_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DNG_INSTALL_DIR_LIB=lib/netgen \
  -DNG_INSTALL_DIR_INCLUDE=include/netgen \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DCPACK_PACKAGING_INSTALL_PREFIX=/usr \
  $CMAKE_ARGS

make -j12
make install
make package
cd ngsolve
if [ "zesty" = "$UBUNTU_VERSION_NAME" ]
then
  echo "build docu"
  make docs
  echo "upload docu"
  rsync -ztrl --del -e ssh \
    docs/html/* \
    gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$CI_PIPELINE_ID/docu/
fi

