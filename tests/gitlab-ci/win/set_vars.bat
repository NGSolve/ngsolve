call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64"
set CMAKE_GENERATOR=Ninja
set NETGEN_BUILD_DIR=%CI_DIR%\build
set CMAKE_INSTALL_PREFIX=%CI_DIR%\install
set NETGENDIR=%CMAKE_INSTALL_PREFIX%\bin
set PATH=%NETGENDIR%;%PATH%
set PYTHONPATH=%CMAKE_INSTALL_PREFIX%\lib\site-packages
set CLCACHE_BASEDIR=%CI_DIR%
set CMAKE_CONFIG=%CMAKE_CONFIG% -DUSE_OCC=ON -DUSE_CCACHE=ON -DUSE_MKL=ON -DMKL_STATIC=ON -DMKL_ROOT="C:/Intel/compilers_and_libraries/windows/mkl"
