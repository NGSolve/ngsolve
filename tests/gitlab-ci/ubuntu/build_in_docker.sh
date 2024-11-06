set -e

# Run ssh-agent (inside the build environment)
eval $(ssh-agent -s)
# Add the SSH key stored in SSH_PRIVATE_KEY variable to the agent store
ssh-add <(echo "$SSH_PRIVATE_KEY")
mkdir -p ~/.ssh
[[ -f /.dockerenv ]] && echo -e "Host *\n\tStrictHostKeyChecking no\n\n" > ~/.ssh/config


if [ "$IMAGE_NAME" == "debug" ]
then
  export CMAKE_CXX_FLAGS="-Og -Wall -Wno-sign-compare -DDebug"
  export CMAKE_BUILD_TYPE="Debug"
  # export CMAKE_ARGS="$CMAKE_ARGS -DFAST_COMPILE=ON"
else
  export CMAKE_BUILD_TYPE="Release"
fi

if [ "$IMAGE_NAME" == "avx" ] || [ "$IMAGE_NAME" == "avx512" ]
then
  export CMAKE_ARGS="$CMAKE_ARGS -DBUILD_OCC=ON"
  export USE_NATIVE_ARCH="ON"
else
  export USE_NATIVE_ARCH="OFF"
fi

if [ "$IMAGE_NAME" == "mpi" ] || [ "$IMAGE_NAME" == "avx" ]
then
  export PYTHONPATH=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real/lib/python3/dist-packages
  export CMAKE_ARGS="$CMAKE_ARGS -DUSE_MPI=ON -DMKL_STATIC=ON -DMKL_SDL=OFF -DUSE_HYPRE=OFF -DUSE_MUMPS=OFF -DMKL_MULTI_THREADED=OFF -DUSE_GUI=OFF -DBUILD_STUB_FILES=OFF"
fi

cd 
cd src/ngsolve
cd external_dependencies
rm -rf netgen
cd ..
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve

mkdir -p /logs/
pip3 freeze > /logs/pip_freeze.log

cmake ../../src/ngsolve \
  -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
  -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
  -DUSE_NATIVE_ARCH=$USE_NATIVE_ARCH \
  -DUSE_CGNS=ON \
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
  $CMAKE_ARGS

make -j12
make install
cd ngsolve

if [ "$IMAGE_NAME" == "avx" ]
then
  ## build and upload docu to server

  # first build ngsxfem
  cd ~/src
  git clone --branch master https://github.com/ngsxfem/ngsxfem.git
  cd ngsxfem
  mkdir -p build
  cd build
  cmake -DBUILD_NGSOLVE=OFF -DCMAKE_INSTALL_PREFIX=/usr ..
  make -j12 install

  # fix links like /edit/some_file.cpp to work with nbsphinx (removing the /edit/ part)
  sed -i 's/\/edit\///g' ~/src/ngsolve/docs/i-tutorials/*/*.ipynb

  cd ~/build/ngsolve/ngsolve
  export NGS_NUM_THREADS=4
  echo "build docu"
  ipython3 profile create --parallel --profile=default
  echo 'c.MPILauncher.mpi_args = ["--allow-run-as-root"]' >> ~/.ipython/profile_default/ipcluster_config.py
  # export SHOW_LOGS=1
  if [[ $SHOW_LOGS ]];
  then
    make docs
  else
    make docs > /logs/build_docs.log 2>&1
  fi
  find ~/src/ngsolve/docs/i-tutorials -name '*.ipynb' -print0 | xargs -0 nbstripout
  cp -r ~/src/ngsolve/docs/i-tutorials docs/html/jupyter-files
  zip -r docs/html/i-tutorials.zip docs/html/jupyter-files
  echo "upload docu"
  rsync -ztrl --del -e ssh \
    --rsync-path="mkdir -p deploy/builds/$CI_PIPELINE_ID/docu/ && rsync" \
    docs/html/* \
    gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$CI_PIPELINE_ID/docu/
fi
