call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64"
set CMAKE_GENERATOR=Ninja
set NETGEN_BUILD_DIR=%CI_DIR%\build
set CMAKE_INSTALL_PREFIX=%CI_DIR%\install
set SRC_DIR=%CI_DIR%\src
set NETGENDIR=%CMAKE_INSTALL_PREFIX%\bin
set PATH=%NETGENDIR%;%PATH%
set PYTHONPATH=%CMAKE_INSTALL_PREFIX%\lib\site-packages
set CLCACHE_BASEDIR=%CI_DIR%
set HOME=%USERPROFILE%
call ccache . -M 20000000000
call ccache . -c
call ccache . -s
