import glob
import os
import sys
import netgen.version
import site

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib
import pkg_resources

setup_requires = []

def install_filter(cmake_manifest):
    return cmake_manifest

def _patched_parse_manifests(self):
    paths = \
        glob.glob(os.path.join(skbuild.cmaker.CMAKE_BUILD_DIR(), "ngsolve", "install_manifest*.txt"))
    try:
        return [self._parse_manifest(path) for path in paths][0]
    except IndexError:
        return []
   
# we are using the ngsolve superbuild (to download and build some dependencies)
# patch the parse_manifests function to point to the actual ngsolve cmake project within the superbuild
skbuild.cmaker.CMaker._parse_manifests = _patched_parse_manifests

version = check_output([sys.executable,'tests/get_python_version_string_from_git.py'], cwd='.').decode('utf-8').strip()

py_install_dir = get_python_lib(1,0,'').replace('\\','/')

netgen_name = netgen.config.NETGEN_PYTHON_PACKAGE_NAME
name = netgen_name.replace("netgen-mesher", "ngsolve") # keep -avx2 suffix

netgen_version = netgen.config.NETGEN_VERSION_PYTHON
netgen_dir = os.path.abspath(os.path.join(netgen.__file__, '../'*(len(py_install_dir.split('/'))+2)))

install_requires = [
        f'{netgen_name} == {netgen_version}',
    ]

_cmake_args = [
    f'-DNETGEN_DIR={netgen_dir}',    
    '-DUSE_SUPERBUILD:BOOL=ON',
    '-DCMAKE_BUILD_TYPE=Release',
    '-DBUILD_FOR_CONDA=ON',
    '-DUSE_UMFPACK=OFF',
    f'-DNGSOLVE_VERSION_PYTHON={version}',
]

if 'NETGEN_CCACHE' in os.environ:
  _cmake_args += ['-DUSE_CCACHE=ON']


if 'darwin' in sys.platform:
    pass
elif 'linux' in sys.platform:
    _cmake_args += [
        '-DUSE_MKL:BOOL=ON',
        f'-DMKL_ROOT:PATH={netgen_dir}',
        f'-DMKL_LIBRARY:PATH={netgen_dir}/lib/libmkl_rt.so.1',
        f'-DMKL_INCLUDE_DIR:PATH={netgen_dir}/include',
    ]
    install_requires.append('mkl == 2021.*')
elif 'win' in sys.platform:
    _cmake_args += [
        '-DUSE_MKL:BOOL=ON',
        f'-DMKL_ROOT:PATH={netgen_dir}',
        f'-DMKL_LIBRARY:PATH={netgen_dir}/Library/lib/mkl_rt.lib',
        f'-DMKL_INCLUDE_DIR:PATH={netgen_dir}/Library/include',
        f'-DNGSOLVE_INSTALL_DIR_TCL:PATH=Scripts',
    ]
    install_requires.append('mkl == 2021.*')

if 'PYDIR' in os.environ:
    _cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}']

setup(
    name=name,
    version=version,
    description="NGSolve",
    author='The NGSolve team',
    license="LGPL2.1",
    packages=['ngsolve'],
    package_dir={'ngsolve': 'python'},
    install_requires=install_requires,
    tests_require=['pytest','scipy','numpy'],
    include_package_data=True,
    cmake_process_manifest_hook=install_filter,
    cmake_args=_cmake_args,
    setup_requires=setup_requires
)
