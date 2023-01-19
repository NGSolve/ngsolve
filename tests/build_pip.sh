#! /bin/bash
set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

cd external_dependencies/netgen && git remote update && cd ../..

rm -rf wheelhouse
mkdir wheelhouse
py=/opt/python/cp39-cp39/bin/python 
export NETGEN_VERSION=`$py tests/get_python_version_string_from_git.py external_dependencies/netgen`
export NETGEN_CCACHE=1

$py tests/fix_auditwheel_policy.py

for pyversion in 38 39 310 311
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2022.* mkl-devel==2022.*
    $PYDIR/pip install netgen-mesher==$NETGEN_VERSION

    rm -rf _skbuild
    $PYDIR/pip wheel .
    rename linux_ manylinux_2_17_x86_64.manylinux2014_ ngsolve*.whl
    mv ngsolve*.whl wheelhouse/
    $PYDIR/pip uninstall -y netgen-mesher

    #$PYDIR/pip install --extra-index-url https://test.pypi.org/simple/ wheelhouse/ngsolve-avx2-*-cp${pyversion}-*.whl
    #$PYDIR/python3 -c 'import ngsolve'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl
