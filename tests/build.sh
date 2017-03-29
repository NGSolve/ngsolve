set -e
cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve -DUSE_NATIVE_ARCH=OFF -DUSE_OCC=ON -DUSE_CCACHE=ON -DUSE_MKL=ON -DUSE_UMFPACK=ON -DINSTALL_PROFILES=ON -DMKL_STATIC=ON -DENABLE_UNIT_TESTS=ON -DCMAKE_BUILD_TYPE=Release -DCPACK_PACKAGING_INSTALL_PREFIX=/opt/netgen $CMAKE_ARGS
make -j12
make install
make package
cd ngsolve
make docs
cd ..

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

if [ "yakkety" = "$UBUNTU_VERSION_NAME" ]
then
  echo "upload docu"
  rsync -ztrl --del -e ssh \
    ngsolve/docs/html/* \
    gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$CI_PIPELINE_ID/docu/
fi

