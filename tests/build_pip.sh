#! /bin/bash
set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

rm -rf wheelhouse
py=/opt/python/cp39-cp39/bin/python 
export NETGEN_VERSION=`$py tests/get_python_version_string_from_git.py external_dependencies/netgen`
export NETGEN_CCACHE=1

$py tests/fix_auditwheel_policy.py

for pyversion in 38 39 310
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2021.* mkl-devel==2021.*
    $PYDIR/pip install -i https://test.pypi.org/simple/ netgen-mesher==$NETGEN_VERSION

    rm -rf _skbuild
    $PYDIR/pip wheel --use-feature=in-tree-build --extra-index-url https://test.pypi.org/simple/ .
    auditwheel repair ngsolve*-cp${pyversion}-*.whl
    rm ngsolve-*.whl
    $PYDIR/pip uninstall -y netgen-mesher


    rm -rf _skbuild
    $PYDIR/pip install -i https://test.pypi.org/simple/ netgen-mesher-avx2==$NETGEN_VERSION
    NETGEN_ARCH=avx2 $PYDIR/pip wheel --use-feature=in-tree-build --extra-index-url https://test.pypi.org/simple/ .
    auditwheel repair ngsolve_avx2*-cp${pyversion}-*.whl
    rm ngsolve_avx2-*.whl

    #$PYDIR/pip install --extra-index-url https://test.pypi.org/simple/ wheelhouse/ngsolve-avx2-*-cp${pyversion}-*.whl
    #$PYDIR/python3 -c 'import ngsolve'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl
