set -e

# Build and upload the macOS wheels for all supported CPython versions using
# cibuildwheel (it downloads the CPython runtimes and runs delocate itself).
rm -rf wheelhouse dist

export PATH=/Applications/CMake.app/Contents/bin:$PATH

export CCACHE_BASEDIR=$(pwd)
export CCACHE_NOHASHDIR=1
export CCACHE_DIR=${CCACHE_DIR:-$HOME/.cache/ccache}
# Bound the cache (all CPython versions in one run) and reset stats; cibuildwheel
# inherits the environment on macOS so these reach the actual build.
if command -v ccache >/dev/null 2>&1; then
  ccache -M ${CCACHE_MAXSIZE:-25G}
  ccache -z
fi

# Netgen version this NGSolve build links against (matching wheel from PyPI).
# Exported so the cibuildwheel before-build step can pin netgen-mesher to it.
export NETGEN_VERSION=`python3 external_dependencies/netgen/tests/utils.py --get-version --dir=./external_dependencies/netgen`
echo "Netgen version: $NETGEN_VERSION"

python3 -m pip install -U pip cibuildwheel twine
CIBW_ENVIRONMENT_PASS_MACOS="NETGEN_VERSION" python3 -m cibuildwheel --platform macos --output-dir wheelhouse
command -v ccache >/dev/null 2>&1 && ccache -s
python3 -m twine upload --skip-existing wheelhouse/*.whl
