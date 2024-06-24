pip3 install --upgrade pybind11-stubgen

git submodule update --init --recursive
IF %errorlevel% NEQ 0 exit /b %errorlevel%
rd /s /q %NETGEN_BUILD_DIR%
mkdir %NETGEN_BUILD_DIR%
IF %errorlevel% NEQ 0 exit /b %errorlevel%
cd %NETGEN_BUILD_DIR%
IF %errorlevel% NEQ 0 exit /b %errorlevel%

cmake %SRC_DIR% ^
        -G"%CMAKE_GENERATOR%" ^
        -DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX% ^
        %CMAKE_CONFIG% ^
        -DUSE_NATIVE_ARCH=%NG_USE_NATIVE_ARCH% ^
        -DUSE_CGNS=ON ^
        -DUSE_OCC=ON ^
        -DUSE_CCACHE=ON ^
        -DUSE_MKL=ON ^
        -DMKL_STATIC=ON ^
        -DMKL_ROOT="C:/Intel/compilers_and_libraries/windows/mkl" ^
        -DUSE_UMFPACK=ON ^
        -DENABLE_UNIT_TESTS=ON ^
        -DCPACK_PACKAGE_NAME=NGSolve ^
        -DCMAKE_BUILD_TYPE=Release

IF %errorlevel% NEQ 0 exit /b %errorlevel%
cmake --build . --target install --config Release
IF %errorlevel% NEQ 0 exit /b %errorlevel%
