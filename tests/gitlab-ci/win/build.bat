git submodule update --init --recursive
IF %errorlevel% NEQ 0 exit /b %errorlevel%
rd /s /q %NETGEN_BUILD_DIR%
mkdir %NETGEN_BUILD_DIR%
IF %errorlevel% NEQ 0 exit /b %errorlevel%
cd %NETGEN_BUILD_DIR%
IF %errorlevel% NEQ 0 exit /b %errorlevel%

cmake %CI_PROJECT_DIR% ^
        -G"%CMAKE_GENERATOR%" ^
        -DCMAKE_INSTALL_PREFIX=%CMAKE_INSTALL_PREFIX% ^
        %CMAKE_CONFIG% ^
        -DUSE_UMFPACK=ON ^
        -DENABLE_UNIT_TESTS=ON ^
        -DCPACK_PACKAGE_NAME=NGSolve ^
        -DCMAKE_BUILD_TYPE=Release

IF %errorlevel% NEQ 0 exit /b %errorlevel%
cmake --build . --target install --config Release
