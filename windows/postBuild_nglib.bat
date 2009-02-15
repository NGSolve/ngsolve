REM *********************************************************************************
REM *** Netgen Library (nglib) Windows Post-Build Script
REM *** Author: Philippose Rajan
REM *** Date: 15/02/2009
REM ***
REM *** Used to perform an "Install" of the generated Dynamic Link Library (DLL)  
REM *** and the corresponding LIB file.
REM ***
REM *** Call from Visual C++ using:
REM *** postBuild_nglib.bat $(ProjectName) $(TargetFileName) $(ConfigurationName) $(ProjectDir)
REM *********************************************************************************
if [%1]==[] goto BuildEventFailed
set PROJ_NAME=%1
set PROJ_EXEC=%2
set BUILD_TYPE=%3
set PROJ_DIR=%4

REM *** Change these Folders if required ***
set NETGEN_TCLSRC=%PROJ_DIR%..\ng
set INSTALL_FOLDER=%PROJ_DIR%%BUILD_TYPE%-lib\

echo POSTBUILD Script for %PROJ_NAME% ........

REM *** Embed the Windows Manifest into the Executable File ***
echo Embedding Manifest into the executable: %PROJ_EXEC% ....
mt.exe -manifest "%PROJ_DIR%%PROJ_NAME%\%BUILD_TYPE%\%PROJ_EXEC%.intermediate.manifest" "-outputresource:%PROJ_DIR%%PROJ_NAME%\%BUILD_TYPE%\%PROJ_EXEC%;2" 
if errorlevel 1 goto BuildEventFailed

REM *** Copy the DLL and LIB Files into the install folder ***
echo Installing required files into %INSTALL_FOLDER% ....
xcopy "%PROJ_DIR%%PROJ_NAME%\%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%" /i /d /y
if errorlevel 1 goto BuildEventFailed
xcopy "%PROJ_DIR%%PROJ_NAME%\%BUILD_TYPE%\%PROJ_NAME%.lib" "%INSTALL_FOLDER%" /i /d /y
if errorlevel 1 goto BuildEventFailed

REM *** Clean up the build directory by deleting the OBJ files ***
REM echo Deleting the %PROJ_NAME% build folder %PROJ_DIR%%PROJ_NAME% ....
REM rmdir %PROJ_DIR%%PROJ_NAME% /s /q

REM *** If there have been no errors so far, we are done ***
goto BuildEventOK
:BuildEventFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Install Manually !!!
exit 1
:BuildEventOK
echo POSTBUILD Script for %PROJ_NAME% completed OK.....!!
