#! /bin/bash
set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance
alias ulimit=echo
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel

curl https://dl.fedoraproject.org/pub/epel/8/Everything/x86_64/Packages/c/ccache-3.7.7-1.el8.x86_64.rpm -o ccache.rpm
dnf -y install ccache.rpm

cd external_dependencies/netgen
git remote update
bash tests/build_pip.sh
cd ../..

rm -rf wheelhouse
mkdir wheelhouse
export NETGEN_VERSION=`/opt/python/cp312-cp312/bin/python external_dependencies/netgen/tests/utils.py --get-version --dir=./external_dependencies/netgen`
export NETGEN_CCACHE=1
echo "Netgen version: $NETGEN_VERSION"

for pyversion in 313 312 311 310 39
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2023.* mkl-devel==2023.*

    # wait until netgen pip package is available
    cd external_dependencies/netgen
    $PYDIR/python3 ./tests/utils.py --wait-pip
    cd ../..

    $PYDIR/python3 ./external_dependencies/netgen/tests/utils.py --check-pip --package ngsolve || continue
    $PYDIR/pip install netgen-mesher==$NETGEN_VERSION

    sed -i 's/set(DLL_EXT ".so")/set(DLL_EXT ".so.2")/' /opt/python/cp${pyversion}-cp${pyversion}/lib/cmake/mkl/MKLConfig.cmake
    rm -rf _skbuild
    $PYDIR/pip wheel --no-clean . || cat /builds/ngsolve/ngsolve/_skbuild/*/cmake-build/dependencies/Stamp/ngsolve/ngsolve-build.log
    rename linux_ manylinux_2_28_x86_64_ ngsolve*.whl
    mv ngsolve*.whl wheelhouse/ || true
    $PYDIR/pip uninstall -y netgen-mesher
    $PYDIR/pip install -U twine
    $PYDIR/twine upload --skip-existing wheelhouse/ngsolve*-cp${pyversion}*manylinux*.whl

    #$PYDIR/pip install --extra-index-url https://test.pypi.org/simple/ wheelhouse/ngsolve-avx2-*-cp${pyversion}-*.whl
    #$PYDIR/python3 -c 'import ngsolve'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done
