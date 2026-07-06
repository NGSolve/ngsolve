def _cmake_to_bool(s):
    return s.upper() not in ['', '0','FALSE','OFF','N','NO','IGNORE','NOTFOUND']

is_python_package    = _cmake_to_bool("@SKBUILD@")

BUILD_STUB_FILES     = _cmake_to_bool("@BUILD_STUB_FILES@")
BUILD_UMFPACK        = _cmake_to_bool("@BUILD_UMFPACK@")
ENABLE_UNIT_TESTS    = _cmake_to_bool("@ENABLE_UNIT_TESTS@")
INSTALL_DEPENDENCIES = _cmake_to_bool("@INSTALL_DEPENDENCIES@")
USE_CCACHE           = _cmake_to_bool("@USE_CCACHE@")
USE_HYPRE            = _cmake_to_bool("@USE_HYPRE@")
USE_LAPACK           = _cmake_to_bool("@USE_LAPACK@")
USE_MKL              = _cmake_to_bool("@USE_MKL@")
USE_MUMPS            = _cmake_to_bool("@USE_MUMPS@")
USE_PARDISO          = _cmake_to_bool("@USE_PARDISO@")
USE_UMFPACK          = _cmake_to_bool("@USE_UMFPACK@")
USE_CUDA             = _cmake_to_bool("@USE_CUDA@")

NETGEN_DIR = "@NETGEN_DIR@"

NGSOLVE_COMPILE_DEFINITIONS         = "@NGSOLVE_COMPILE_DEFINITIONS@"
NGSOLVE_COMPILE_DEFINITIONS_PRIVATE = "@NGSOLVE_COMPILE_DEFINITIONS_PRIVATE@"
NGSOLVE_COMPILE_INCLUDE_DIRS        = "@NGSOLVE_COMPILE_INCLUDE_DIRS@"
NGSOLVE_COMPILE_OPTIONS             = "@NGSOLVE_COMPILE_OPTIONS@"

NGSOLVE_INSTALL_DIR_PYTHON   = "@NGSOLVE_INSTALL_DIR_PYTHON@"
NGSOLVE_INSTALL_DIR_BIN      = "@NGSOLVE_INSTALL_DIR_BIN@"
NGSOLVE_INSTALL_DIR_LIB      = "@NGSOLVE_INSTALL_DIR_LIB@"
NGSOLVE_INSTALL_DIR_INCLUDE  = "@NGSOLVE_INSTALL_DIR_INCLUDE@"
NGSOLVE_INSTALL_DIR_CMAKE    = "@NGSOLVE_INSTALL_DIR_CMAKE@"
NGSOLVE_INSTALL_DIR_RES      = "@NGSOLVE_INSTALL_DIR_RES@"

NGSOLVE_VERSION = "@NGSOLVE_VERSION@"
NGSOLVE_VERSION_GIT = "@git_version_string@"
NGSOLVE_VERSION_PYTHON = "@NGSOLVE_VERSION_PYTHON@"

NGSOLVE_VERSION_MAJOR = "@NGSOLVE_VERSION_MAJOR@"
NGSOLVE_VERSION_MINOR = "@NGSOLVE_VERSION_MINOR@"
NGSOLVE_VERSION_TWEAK = "@NGSOLVE_VERSION_TWEAK@"
NGSOLVE_VERSION_PATCH = "@NGSOLVE_VERSION_PATCH@"
NGSOLVE_VERSION_HASH = "@NGSOLVE_VERSION_HASH@"

CMAKE_CXX_COMPILER           = "@CMAKE_CXX_COMPILER@"
CMAKE_CUDA_COMPILER          = "@CMAKE_CUDA_COMPILER@"
CMAKE_C_COMPILER             = "@CMAKE_C_COMPILER@"
CMAKE_LINKER                 = "@CMAKE_LINKER@"
CMAKE_INSTALL_PREFIX         = "@CMAKE_INSTALL_PREFIX@"
CMAKE_CXX_COMPILER_LAUNCHER  = "@CMAKE_CXX_COMPILER_LAUNCHER@"

version = NGSOLVE_VERSION_GIT

MKL_LINK = "@MKL_LINK@"

def get_cmake_dir():
    import os.path as p
    d_python = p.dirname(p.dirname(p.dirname(__file__)))
    py_to_cmake = p.relpath(
            NGSOLVE_INSTALL_DIR_CMAKE,
            NGSOLVE_INSTALL_DIR_PYTHON
            )
    return p.normpath(p.join(d_python,py_to_cmake))


def _resolve_install_dir(install_dir, python_dir, anchor_file):
    import os.path as p
    d_python = p.dirname(p.dirname(p.dirname(anchor_file)))
    return p.normpath(p.join(d_python, p.relpath(install_dir, python_dir)))


def get_include_dir():
    "Absolute path of the installed NGSolve headers."
    return _resolve_install_dir(
        NGSOLVE_INSTALL_DIR_INCLUDE, NGSOLVE_INSTALL_DIR_PYTHON, __file__)


def get_library_dir():
    "Absolute path of the installed NGSolve libraries."
    return _resolve_install_dir(
        NGSOLVE_INSTALL_DIR_LIB, NGSOLVE_INSTALL_DIR_PYTHON, __file__)


def get_netgen_include_dir():
    "Absolute path of the installed Netgen headers."
    import netgen.config as ngc
    return _resolve_install_dir(
        ngc.NG_INSTALL_DIR_INCLUDE, ngc.NG_INSTALL_DIR_PYTHON, ngc.__file__)


def get_netgen_library_dir():
    "Absolute path of the installed Netgen libraries."
    import netgen.config as ngc
    return _resolve_install_dir(
        ngc.NG_INSTALL_DIR_LIB, ngc.NG_INSTALL_DIR_PYTHON, ngc.__file__)


def get_include_dirs():
    "All include directories needed to compile against NGSolve (NGSolve and Netgen)."
    import os.path as p
    netgen_inc = get_netgen_include_dir()
    dirs = [get_include_dir(), netgen_inc, p.join(netgen_inc, "include")]
    unique = []
    for d in dirs:
        if d not in unique:
            unique.append(d)
    return unique


def get_library_dirs():
    "All library directories needed to link against NGSolve (NGSolve and Netgen)."
    dirs = [get_library_dir(), get_netgen_library_dir()]
    unique = []
    for d in dirs:
        if d not in unique:
            unique.append(d)
    return unique


def get_compile_include_args():
    "Compiler include flags (NGSolve + Netgen headers), formatted for the current platform."
    import sys
    if sys.platform == "win32":
        q = chr(34)
        return " ".join("/I" + q + d + q for d in get_include_dirs())
    return " ".join("-I" + d for d in get_include_dirs())


def get_link_lib_args():
    "Linker flags to link a shared library against NGSolve and Netgen, for the current platform."
    import sys
    dirs = get_library_dirs()
    if sys.platform == "win32":
        q = chr(34)
        args = ["/LIBPATH:" + q + d + q for d in dirs]
        args += ["nglib.lib", "ngcore.lib", "libngsolve.lib"]
        return " ".join(args)
    if sys.platform == "darwin":
        args = ["-L" + d for d in dirs]
        args.append("-undefined dynamic_lookup")
        return " ".join(args)
    args = []
    for d in dirs:
        args.append("-L" + d)
        args.append("-Wl,-rpath," + d)
    return " ".join(args)
