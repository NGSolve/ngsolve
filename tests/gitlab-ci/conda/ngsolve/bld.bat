mkdir build
cd build

set PYTHONPATH=%PREFIX%\Lib\site-packages

cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ^
       -DUSE_SUPERBUILD=OFF ^
       -DCMAKE_SYSTEM_PREFIX_PATH=%LIBRARY_PREFIX% ^
       -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
       -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
       -DUSE_MKL=ON ^
       -DMKL_SDL=ON ^
       -DMKL_INCLUDE_DIR=%LIBRARY_PREFIX%/include ^
       -DMKL_MINIMAL_LIBRARY=%LIBRARY_PREFIX%/lib/libmkl_rt.lib ^
       -DMKL_LIBRARY=%LIBRARY_PREFIX%/lib/libmkl_rt.lib ^
       -DMKL_ROOT=%LIBRARY_PREFIX% ^
       -DUSE_UMFPACK=OFF ^
       %SRC_DIR%
if errorlevel 1 exit 1

REM Build
nmake
if errorlevel 1 exit 1

REM Install
nmake install
if errorlevel 1 exit 1
