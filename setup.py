import glob
import os
import sys
import netgen.version

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib
from urllib.request import urlopen
import json

setup_requires = ['pybind11-stubgen==2.5']

def install_filter(cmake_manifest):
    return cmake_manifest

def _patched_parse_manifests(self):
    paths = \
        glob.glob(os.path.join(skbuild.cmaker.CMAKE_BUILD_DIR(), "ngsolve", "install_manifest*.txt"))
    try:
        return [self._parse_manifest(path) for path in paths][0]
    except IndexError:
        return []

def get_openblas_requirements():
    if 'darwin' in sys.platform:
        return []

    if 'linux' in sys.platform or 'win' in sys.platform:
        import importlib.metadata
        import ngsolve_openblas
        version = importlib.metadata.version('ngsolve-openblas')
        return [f"ngsolve-openblas=={version}"]

    return []

# we are using the ngsolve superbuild (to download and build some dependencies)
# patch the parse_manifests function to point to the actual ngsolve cmake project within the superbuild
skbuild.cmaker.CMaker._parse_manifests = _patched_parse_manifests

utils_command = [sys.executable, os.path.join('external_dependencies', 'netgen', 'tests', 'utils.py')]
git_version = check_output([*utils_command, '--get-git-version']).decode('utf-8').strip()
version = check_output([*utils_command, '--get-version']).decode('utf-8').strip()

# check if release alread exists on pypi
try:
    python_version = f"cp{sys.version_info.major}{sys.version_info.minor}"
    platform = sys.platform
    if platform == 'darwin':
        platform = 'macosx'
    url = f"https://pypi.org/pypi/ngsolve/{version}/json"
    data = json.loads(urlopen(url).read())
    for url in data['urls']:
        if platform in url['filename'] and url['python_version'] == python_version:
            print("version already exists on pypi, skip build")
            sys.exit()
except:
    pass


py_install_dir = get_python_lib(1,0,'').replace('\\','/')

netgen_name = netgen.config.NETGEN_PYTHON_PACKAGE_NAME
name = netgen_name.replace("netgen-mesher", "ngsolve") # keep -avx2 suffix

netgen_version = netgen.config.NETGEN_VERSION_PYTHON
root_dir = os.path.abspath(os.path.join(netgen.__file__, '../'*(len(py_install_dir.split('/'))+2)))

if netgen.config.NG_INSTALL_DIR_CMAKE.startswith('netgen'):
    netgen_dir = os.path.abspath(os.path.join(os.path.dirname(netgen.__file__), 'cmake'))
else:
    netgen_dir = root_dir

install_requires = [
        f'{netgen_name} == {netgen_version}',
    ] + get_openblas_requirements()

_cmake_args = [
    f'-DNETGEN_DIR={netgen_dir}',
    '-DUSE_SUPERBUILD:BOOL=ON',
    '-DCMAKE_BUILD_TYPE=Release',
    '-DBUILD_FOR_CONDA=ON',
    '-DUSE_UMFPACK=ON',
    '-DBUILD_UMFPACK=ON',
    f'-DNGSOLVE_VERSION_PYTHON={version}',
]

if get_openblas_requirements():
    import ngsolve_openblas
    libs = ";".join(ngsolve_openblas.get_libraries())
    libs = libs.replace('\\', '/')
    _cmake_args += [f'-DLAPACK_LIBRARIES={libs}']


if 'NETGEN_CCACHE' in os.environ:
  _cmake_args += ['-DUSE_CCACHE=ON']

packages=['netgen', 'ngsolve']

if 'darwin' in sys.platform:
    _cmake_args += [
        '-DBUILD_STUB_FILES=ON',
    ]
elif 'linux' in sys.platform:
    _cmake_args += [
        '-DUSE_CUDA=OFF',
        '-DCMAKE_CUDA_ARCHITECTURES=all',
        '-DBUILD_STUB_FILES=ON',
    ]
    packages = []
elif 'win' in sys.platform:
    _cmake_args += [
        f'-DNGSOLVE_INSTALL_DIR_TCL:PATH=Scripts',
        '-DBUILD_STUB_FILES=OFF',
    ]

if 'PYDIR' in os.environ:
    _cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}']

setup(
    name=name,
    version=version,
    description="NGSolve",
    author='The NGSolve team',
    license="LGPL2.1",
    packages=packages,
    # package_dir={'ngsolve': 'python'},
    install_requires=install_requires,
    tests_require=['pytest','scipy','numpy'],
    # include_package_data=True,
    cmake_process_manifest_hook=install_filter,
    cmake_args=_cmake_args,
    setup_requires=setup_requires
)
