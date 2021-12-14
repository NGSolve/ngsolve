mkdir build
cd build

set PYTHONPATH=%PREFIX%\Lib\site-packages

cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ^
       -DCMAKE_SYSTEM_PREFIX_PATH=%LIBRARY_PREFIX% ^
       -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
       -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
       -DUSE_SUPERBUILD=OFF ^
       -DUSE_NATIVE_ARCH=OFF ^
       -DUSE_GUI=ON ^
       -DUSE_JPEG=ON ^
       -DUSE_OCC=ON ^
       -DOCC_INCLUDE_DIR=%LIBRARY_PREFIX%/include/opencascade ^
       -DUSE_MPEG=OFF ^
       -DNG_INSTALL_DIR_LIB=Library/lib/netgen ^
       -DNG_INSTALL_DIR_BIN=Library/bin ^
       -DNG_INSTALL_DIR_RES=Library/lib/share ^
       -DNG_INSTALL_DIR_CMAKE=Library/cmake ^
       -DNG_INSTALL_DIR_INCLUDE=Library/include/netgen ^
       -DTCL_LIBRARY:FILEPATH=%LIBRARY_PREFIX%/lib/tcl86t.lib ^
       -DTK_LIBRARY:FILEPATH=%LIBRARY_PREFIX%/lib/tk86t.lib ^
       -DTK_INCLUDE_PATH=%SRC_DIR%/tk/generic ^
       -DTCL_INCLUDE_PATH=%SRC_DIR%/tcl/generic ^
       -DPYBIND11_INCLUDE_DIR:FILEPATH=%LIBRARY_PREFIX%/include ^
       -DPYTHON_INCLUDE_DIRS:FILEPATH=%LIBRARY_PREFIX%/include ^
       -DNG_INSTALL_PYBIND=OFF ^
       %SRC_DIR%
if errorlevel 1 exit 1

REM Build
nmake
if errorlevel 1 exit 1

REM Install
nmake install
if errorlevel 1 exit 1
