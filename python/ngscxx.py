from shutil import which
import glob
import pathlib
import subprocess
import sys
import sysconfig

import netgen.config as ngc
import ngsolve.config as ngsc

def find_compiler():
    compiler = which(ngsc.CMAKE_CXX_COMPILER) or which(ngsc.CMAKE_CXX_COMPILER).split('/')[-1] or which('g++') or which('c++') or which('clang++')
    if compiler is None:
        raise RuntimeError("No compiler found")
    ccache = which('ccache')
    if ccache:
        return [ccache, compiler]
    else:
        return [compiler]

def find_includes():
    for prefix in [ngc.CMAKE_INSTALL_PREFIX]:
        files = sorted(pathlib.Path(prefix).glob('include/nglib.h'))
        if len(files):
            return "-I"+str(files[0].parent)

if __name__ == '__main__':
    compiler = find_compiler()
    flags = []
    flags.append(ngc.ngcore_compile_options)
    for d in ngc.ngcore_compile_definitions.split(';') + ngsc.NGSOLVE_COMPILE_DEFINITIONS.split(';'):
        flags.append(f"-D{d}")

    for o in ngsc.NGSOLVE_COMPILE_OPTIONS.split(';'):
        if o.startswith("$<$<COMPILE_LANGUAGE:CXX>:"):
            o = o.split(':')[2].split('>')[0]
        flags.append(o)

    flags.append(find_includes())
    flags.append(f"-I{sysconfig.get_path('include')}")

    subprocess.run(compiler + flags + sys.argv[1:])
