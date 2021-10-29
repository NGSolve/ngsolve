set -e
rm -rf _skbuild dist ../venv_ngs

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_CCACHE=1

export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin
export NETGEN_VERSION=`$PYDIR/python3 tests/get_python_version_string_from_git.py external_dependencies/netgen`

$PYDIR/python3 --version
$PYDIR/python3 -m venv ../venv_ngs

source ../venv_ngs/bin/activate
pip install numpy twine scikit-build wheel

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
pip install -i https://test.pypi.org/simple/ netgen-mesher==$NETGEN_VERSION
python setup.py bdist_wheel --plat-name macosx-10.14-universal2 -j10
python -m twine upload --repository testpypi dist/*.whl
