set -e
rm -rf _skbuild dist ../venv_ngs

export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin

cd external_dependencies/netgen
git remote update
bash tests/build_pip_mac.sh $1
$PYDIR/python3 tests/utils.py --wait-pip
cd ../..

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_CCACHE=1

export NETGEN_VERSION=`$PYDIR/python3 external_dependencies/netgen/tests/utils.py --get-version --dir=./external_dependencies/netgen`
echo "Netgen version: $NETGEN_VERSION"

$PYDIR/python3 --version
$PYDIR/python3 external_dependencies/netgen/tests/utils.py --check-pip --package ngsolve || exit 0
$PYDIR/python3 -m venv ../venv_ngs

source ../venv_ngs/bin/activate
$PYDIR/pip3 install numpy twine scikit-build wheel

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
$PYDIR/pip3 install netgen-mesher==$NETGEN_VERSION
$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -j10
$PYDIR/python3 -m twine upload dist/*.whl || true
