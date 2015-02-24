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
   
REM *** Start the Installation procedure ***
echo POSTBUILD Script for %PROJ_NAME% ........


REM *** Copy the primary libmyngsolve DLL into the Install Folder ***
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


REM *** Copy the NgSolve pde tutorials to the install folder ***
echo Copying pde tutorials
C:\Windows\System32\xcopy "%PROJ_DIR%..\*.pde" "%INSTALL_FOLDER%\my_little_ngsolve\*.pde" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\*.vol" "%INSTALL_FOLDER%\my_little_ngsolve\*.vol" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\*.in2d" "%INSTALL_FOLDER%\my_little_ngsolve\*.in2d" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\*.vol.*" "%INSTALL_FOLDER%\my_little_ngsolve\*.vol.*" /i /d /y

C:\Windows\System32\xcopy "%PROJ_DIR%..\mixed\*.pde" "%INSTALL_FOLDER%\my_little_ngsolve\mixed\*.pde" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\mixed\*.vol" "%INSTALL_FOLDER%\my_little_ngsolve\mixed\*.vol" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\mixed\*.in2d" "%INSTALL_FOLDER%\my_little_ngsolve\mixed\*.in2d" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\mixed\*.vol.*" "%INSTALL_FOLDER%\my_little_ngsolve\mixed\*.vol.*" /i /d /y

C:\Windows\System32\xcopy "%PROJ_DIR%..\precond\*.pde" "%INSTALL_FOLDER%\my_little_ngsolve\precond\*.pde" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\precond\*.vol" "%INSTALL_FOLDER%\my_little_ngsolve\precond\*.vol" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\precond\*.in2d" "%INSTALL_FOLDER%\my_little_ngsolve\precond\*.in2d" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\precond\*.vol.*" "%INSTALL_FOLDER%\my_little_ngsolve\precond\*.vol.*" /i /d /y

C:\Windows\System32\xcopy "%PROJ_DIR%..\tcl\*.pde" "%INSTALL_FOLDER%\my_little_ngsolve\tcl\*.pde" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\tcl\*.vol" "%INSTALL_FOLDER%\my_little_ngsolve\tcl\*.vol" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\tcl\*.in2d" "%INSTALL_FOLDER%\my_little_ngsolve\tcl\*.in2d" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\tcl\*.vol.*" "%INSTALL_FOLDER%\my_little_ngsolve\tcl\*.vol.*" /i /d /y

C:\Windows\System32\xcopy "%PROJ_DIR%..\HDG\*.pde" "%INSTALL_FOLDER%\my_little_ngsolve\HDG\*.pde" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\HDG\*.vol" "%INSTALL_FOLDER%\my_little_ngsolve\HDG\*.vol" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\HDG\*.in2d" "%INSTALL_FOLDER%\my_little_ngsolve\HDG\*.in2d" /i /d /y
C:\Windows\System32\xcopy "%PROJ_DIR%..\HDG\*.vol.*" "%INSTALL_FOLDER%\my_little_ngsolve\HDG\*.vol.*" /i /d /y

echo %PROJ_NAME%
if /i "%PROJ_NAME%" == "libtcldemo" (
C:\Windows\System32\xcopy "%PROJ_DIR%..\tcl\*.tcl" "%INSTALL_FOLDER%\bin\*.tcl" /i /d /y
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
:ExecInstallFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying the %PROJ_NAME% executable into install folder!!!
exit 1
:LibInstallFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying %LIB_NAME%.lib or %LIB_NAME%.h into install folder!!!
exit 1

:BuildEventOK
echo POSTBUILD Script for %PROJ_NAME% completed OK.....!!

