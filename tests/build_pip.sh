#! /bin/bash
set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance
sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-*
sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

cd external_dependencies/netgen
git remote update
bash tests/build_pip.sh
cd ../..

rm -rf wheelhouse
mkdir wheelhouse
py=/opt/python/cp39-cp39/bin/python 
export NETGEN_VERSION=`$py tests/get_python_version_string_from_git.py external_dependencies/netgen`
export NETGEN_CCACHE=1

$py tests/fix_auditwheel_policy.py

for pyversion in 312 311 310 39 38
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2023.* mkl-devel==2023.*
    $PYDIR/pip install netgen-mesher==$NETGEN_VERSION

    sed -i 's/set(DLL_EXT ".so")/set(DLL_EXT ".so.2")/' /opt/python/cp${pyversion}-cp${pyversion}/lib/cmake/mkl/MKLConfig.cmake
    rm -rf _skbuild
    $PYDIR/pip wheel --no-clean . || cat /builds/ngsolve/ngsolve/_skbuild/*/cmake-build/dependencies/Stamp/ngsolve/ngsolve-build.log
    rename linux_ manylinux_2_17_x86_64.manylinux2014_ ngsolve*.whl
    mv ngsolve*.whl wheelhouse/ || true
    $PYDIR/pip uninstall -y netgen-mesher
    $PYDIR/pip install -U twine
    $PYDIR/twine upload --skip-existing wheelhouse/ngsolve*-cp${pyversion}*manylinux*.whl

    #$PYDIR/pip install --extra-index-url https://test.pypi.org/simple/ wheelhouse/ngsolve-avx2-*-cp${pyversion}-*.whl
    #$PYDIR/python3 -c 'import ngsolve'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done
