REM *********************************************************************************
REM *** NGSolve Windows Post-Build Script
REM *** Author: Philippose Rajan
REM *** Date: 02/03/2010
REM ***
REM *** Used to perform an "Install" of the generated DLLs and Libraries 
REM *** along with the required *.tcl files
REM ***
REM *** Call from Visual C++ using:
REM *** postBuild_ngsolve.bat $(ProjectName) $(TargetFileName) $(ConfigurationName)  $(PlatformName) $(ProjectDir) <lib filename>
REM *********************************************************************************
if [%1]==[] goto InputParamsFailed
set PROJ_NAME=%~1
set PROJ_EXEC=%~2
set BUILD_TYPE=%~3
set BUILD_ARCH=%~4
set PROJ_DIR=%~5
set LIB_NAME=%~6


REM *** Change these Folders if required ***
REM Check if the environment variable NETGENDIR exists, 
REM and use it as the installation folder
if defined NETGENDIR (
   echo Environment variable NETGENDIR found: %NETGENDIR%
   set INSTALL_FOLDER=%NETGENDIR%\..
) else (
   echo Environment variable NETGENDIR not found.... using default location!!!
   set INSTALL_FOLDER=%PROJ_DIR%..\..\%PROJ_NAME%-inst_%BUILD_ARCH%
)
   
set NGSOLVE_TCLSRC=%PROJ_DIR%..
set NGSOLVE_INCROOT=%PROJ_DIR%..

REM *** Start the Installation procedure ***
echo POSTBUILD Script for %PROJ_NAME% ........

REM *** Copy the core TCL files into the Install Folder ***
echo Installing core TCL files into %INSTALL_FOLDER%\bin ....
C:\Windows\System32\xcopy  "%NGSOLVE_TCLSRC%\*.tcl" "%INSTALL_FOLDER%\bin\" /i /d /y
if errorlevel 1 goto CoreTCLFailed
echo Installing core TCL Files: Completed OK!!

REM *** Copy the primary NGSolve DLL into the Install Folder ***
echo Installing %PROJ_EXEC% into %INSTALL_FOLDER%\bin ....
if /i "%BUILD_ARCH%" == "win32" (
   C:\Windows\System32\xcopy "%PROJ_DIR%%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%\bin\" /i /d /y
   if errorlevel 1 goto ExecInstallFailed
)
if /i "%BUILD_ARCH%" == "x64" (
   C:\Windows\System32\xcopy "%PROJ_DIR%%BUILD_ARCH%\%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%\bin\" /i /d /y
   if errorlevel 1 goto ExecInstallFailed
)   
echo Installing %PROJ_EXEC%: Completed OK!!

REM *** Copy the NGSolve library and include files into the Install Folder ***
echo Installing %LIB_NAME%.lib into %INSTALL_FOLDER%\lib ....
if /i "%BUILD_ARCH%" == "win32" (
   C:\Windows\System32\xcopy "%PROJ_DIR%%BUILD_TYPE%\%LIB_NAME%.lib" "%INSTALL_FOLDER%\lib\" /i /d /y
   if errorlevel 1 goto LibInstallFailed
)
if /i "%BUILD_ARCH%" == "x64" (
   C:\Windows\System32\xcopy "%PROJ_DIR%%BUILD_ARCH%\%BUILD_TYPE%\%LIB_NAME%.lib" "%INSTALL_FOLDER%\lib\" /i /d /y
   if errorlevel 1 goto LibInstallFailed
)   
echo Installing %LIB_NAME%.lib: Completed OK!!

echo Installing NgSolve development header files into %INSTALL_FOLDER%\include ....
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\include\*.hpp" "%INSTALL_FOLDER%\include\include\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\solve\*.hpp" "%INSTALL_FOLDER%\include\solve\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\parallel\*.hpp" "%INSTALL_FOLDER%\include\parallel\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\ngstd\*.hpp" "%INSTALL_FOLDER%\include\ngstd\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\multigrid\*.hpp" "%INSTALL_FOLDER%\include\multigrid\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\linalg\*.hpp" "%INSTALL_FOLDER%\include\linalg\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\fem\*.hpp" "%INSTALL_FOLDER%\include\fem\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\comp\*.hpp" "%INSTALL_FOLDER%\include\comp\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\basiclinalg\*.hpp" "%INSTALL_FOLDER%\include\basiclinalg\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\basiclinalg\*.h" "%INSTALL_FOLDER%\include\basiclinalg\" /i /d /y
C:\Windows\System32\xcopy "%NGSOLVE_INCROOT%\ngstd\*.py" "%INSTALL_FOLDER%\bin\" /i /d /y


REM *** Copy the NgSolve pde tutorials to the install folder ***
echo Copying pde tutorials
C:\Windows\System32\xcopy "%PROJ_DIR%..\pde_tutorial" "%INSTALL_FOLDER%\pde_tutorial\" /i /d /y /e

REM *** Copy the Lapack libraries to the bin folder ***
if /i "%BUILD_TYPE%"=="Release_Lapack" (
if /i "%BUILD_ARCH%"=="x64"	(
		echo Copying Lapack x64 dlls and gfortran runtime libraries
		C:\Windows\System32\xcopy "%PROJ_DIR%..\..\ext_libs\OpenBlas\lib\x64\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y
)
if /i "%BUILD_ARCH%"=="win32" (
		echo Copying Lapack win32 dlls and gfortran runtimes libraries
		C:\Windows\System32\xcopy "%PROJ_DIR%..\..\ext_libs\OpenBlas\lib\x86\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y
)	
)



REM *** Done with the installation routine ***

REM *** Clean up the build directory by deleting the OBJ files ***
REM echo Deleting the %PROJ_NAME% build folder %PROJ_DIR%%PROJ_NAME% ....
REM rmdir %PROJ_DIR%%BUILD_TYPE% /s /q

REM *** If there have been no errors so far, we are done ***
goto BuildEventOK

REM *** Error Messages for each stage of the post build process ***
:InputParamsFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Invalid number of input parameters!!!
exit 1
:ManifestFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Manifest not successfully embedded!!!
exit 1
:CoreTCLFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying core TCL Files into install folder!!!
exit 1
:ExecInstallFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying the %PROJ_NAME% executable into install folder!!!
exit 1
:LibInstallFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying %LIB_NAME%.lib or %LIB_NAME%.h into install folder!!!
exit 1

:BuildEventOK
echo POSTBUILD Script for %PROJ_NAME% completed OK.....!!

