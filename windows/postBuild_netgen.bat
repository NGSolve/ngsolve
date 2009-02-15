REM *********************************************************************************
REM *** Netgen Windows Post-Build Script
REM *** Author: Philippose Rajan
REM *** Date: 31/01/2009
REM ***
REM *** Used to perform an "Install" of the generated executable 
REM *** along with the required *.tcl files
REM ***
REM *** Call from Visual C++ using:
REM *** postBuild_netgen.bat $(ProjectName) $(TargetFileName) $(ConfigurationName) $(ProjectDir)
REM *********************************************************************************
if [%1]==[] goto BuildEventFailed
set PROJ_NAME=%1
set PROJ_EXEC=%2
set BUILD_TYPE=%3
set PROJ_DIR=%4

REM *** Change these Folders if required ***
set NETGEN_TCLSRC=%PROJ_DIR%..\ng
set INSTALL_FOLDER=%PROJ_DIR%%BUILD_TYPE%-bin

echo POSTBUILD Script for %PROJ_NAME% ........

REM *** Embed the Windows Manifest into the Executable File ***
echo Embedding Manifest into the executable: %PROJ_EXEC% ....
mt.exe -manifest "%BUILD_TYPE%\%PROJ_EXEC%.intermediate.manifest" "-outputresource:%BUILD_TYPE%\%PROJ_EXEC%;1" 
if errorlevel 1 goto BuildEventFailed

REM *** Copy the TCL files and the executable to the Install Folder ***
echo Installing required files into %INSTALL_FOLDER% ....
xcopy  "%NETGEN_TCLSRC%\*.tcl" "%INSTALL_FOLDER%" /i /d /y
if errorlevel 1 goto BuildEventFailed

xcopy "%PROJ_DIR%%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%" /i /d /y
if errorlevel 1 goto BuildEventFailed

REM *** Clean up the build directory by deleting the OBJ files ***
REM echo Deleting the %PROJ_NAME% build folder %PROJ_DIR%%PROJ_NAME% ....
REM rmdir %PROJ_DIR%%BUILD_TYPE% /s /q

REM *** If there have been no errors so far, we are done ***
goto BuildEventOK
:BuildEventFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Install Manually !!!
exit 1
:BuildEventOK
echo POSTBUILD Script for %PROJ_NAME% completed OK.....!!
