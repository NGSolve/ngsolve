cd %NETGEN_BUILD_DIR%
IF %errorlevel% NEQ 0 exit /b %errorlevel%
cmake --build . --target package --config Release
IF %errorlevel% NEQ 0 exit /b %errorlevel%
set HOME=/cygdrive/c/tools
rsync -ztrl --del -e "ssh -o StrictHostKeyChecking=no" ^
      --rsync-path="mkdir -p deploy/builds/%CI_PIPELINE_ID%/windows && rsync" ^
      *.msi ^
      gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/%CI_PIPELINE_ID%/windows/

IF %errorlevel% NEQ 0 exit /b %errorlevel%
