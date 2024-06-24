call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64"
set CMAKE_GENERATOR=Ninja
set NETGEN_BUILD_DIR=%CI_PROJECT_DIR%\build
set CMAKE_INSTALL_PREFIX=%CI_PROJECT_DIR%\install
set SRC_DIR=%CI_PROJECT_DIR%
set NETGENDIR=%CMAKE_INSTALL_PREFIX%\bin
set PATH=%NETGENDIR%;C:\python312;C:\python312\bin;C:\python312\Scripts;C:\tools\;%PATH%
set PYTHONPATH=%CMAKE_INSTALL_PREFIX%\lib\site-packages
REM set CCACHE_BASEDIR=%CI_DIR%
set HOME=%USERPROFILE%
call ccache -M 20G
call ccache -s
