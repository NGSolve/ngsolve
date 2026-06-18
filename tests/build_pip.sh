#! /bin/bash
set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance

# This script runs *inside* the manylinux build image (see pip_linux in
# tests/gitlab-ci/pip.yml)

export NETGEN_VERSION=`/opt/python/cp312-cp312/bin/python external_dependencies/netgen/tests/utils.py --get-version --dir=./external_dependencies/netgen`
echo "Netgen version: $NETGEN_VERSION"

rm -rf wheelhouse dist
mkdir -p wheelhouse
export NETGEN_ARCH=avx2

export CCACHE_BASEDIR=$(pwd)
export CCACHE_NOHASHDIR=1
export CCACHE_DIR=${CCACHE_DIR:-$HOME/.ccache}

for pyversion in 314 313 312 311 310
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pip requests packaging
    # skip this version if it is already published on PyPI
    $PYDIR/python3 ./tests/utils.py --check-pip --package ngsolve || continue

    $PYDIR/pip install -U build twine ninja cmake "scikit-build-core>=0.10" pybind11-stubgen
    $PYDIR/pip install -U netgen-occt==7.8.1 netgen-occt-devel==7.8.1 netgen-mesher==$NETGEN_VERSION "ngsolve-openblas==$OPENBLAS_VERSION_PIP"

    # Point CMake at the Python prefix for Tcl/Tk discovery and pin the python root
    PYPREFIX=$($PYDIR/python3 -c 'import sys; print(sys.prefix)')
    export CMAKE_ARGS="-DCMAKE_PREFIX_PATH=${PYPREFIX} -DPython3_ROOT_DIR=${PYPREFIX}"

    rm -rf dist
    $PYDIR/python3 -m build --wheel --no-isolation --outdir dist .

    ccache -s

    # relabel linux_x86_64 -> manylinux (see comment above)
    rename linux_ manylinux_2_28_ dist/ngsolve*-cp${pyversion}-*.whl
    mv dist/ngsolve*-cp${pyversion}-*.whl wheelhouse/

    $PYDIR/pip install wheelhouse/ngsolve*-cp${pyversion}-*manylinux*.whl
    $PYDIR/python3 -c 'import ngsolve; print(ngsolve.__version__)'
    $PYDIR/twine upload --skip-existing wheelhouse/ngsolve*-cp${pyversion}*manylinux*.whl
done
