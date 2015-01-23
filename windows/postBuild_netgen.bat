REM *********************************************************************************
REM *** Netgen Windows Post-Build Script
REM *** Author: Philippose Rajan
REM *** Date: 31/01/2009
REM ***
REM *** Used to perform an "Install" of the generated executable 
REM *** along with the required *.tcl files
REM ***
REM *** Call from Visual C++ using:
REM *** postBuild_netgen.bat $(ProjectName) $(TargetFileName) $(ConfigurationName) $(ProjectDir) <lib filename>
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
set W_WO_OCC=%BUILD_TYPE:~-4,3%
if defined NETGENDIR (
   echo Environment variable NETGENDIR found: %NETGENDIR%
   set INSTALL_FOLDER=%NETGENDIR%\..
) else (
   echo Environment variable NETGENDIR not found.... using default location!!!
   if /i "%W_WO_OCC%" == "OCC" (
      set INSTALL_FOLDER=%PROJ_DIR%..\..\%PROJ_NAME%-instOCC_%BUILD_ARCH%
   ) else (
      set INSTALL_FOLDER=%PROJ_DIR%..\..\%PROJ_NAME%-instNoOCC_%BUILD_ARCH%
   )   
)
   
set NETGEN_TCLSRC=%PROJ_DIR%..\ng
set NETGEN_LIBINC=%PROJ_DIR%..\libsrc\include
set NETGEN_NGSINC=%PROJ_DIR%..\libsrc

REM *** Start the Installation procedure ***
echo POSTBUILD Script for %PROJ_NAME% ........

REM *** Embed the Windows Manifest into the Executable File ***
REM echo Embedding Manifest into the executable: %PROJ_EXEC% ....
REM mt.exe -manifest "%BUILD_TYPE%\%PROJ_EXEC%.intermediate.manifest" "-outputresource:%BUILD_TYPE%\%PROJ_EXEC%;#1" 
REM if errorlevel 1 goto ManifestFailed
REM echo Embedding Manifest into the executable: Completed OK!!

REM *** Copy the core TCL files into the Install Folder ***
echo Installing core TCL files into %INSTALL_FOLDER%\bin ....
xcopy  "%NETGEN_TCLSRC%\*.tcl" "%INSTALL_FOLDER%\bin\" /i /d /y
xcopy  "%NETGEN_TCLSRC%\*.ocf" "%INSTALL_FOLDER%\bin\" /i /d /y
if errorlevel 1 goto CoreTCLFailed
echo Installing core TCL Files: Completed OK!!

REM *** Copy any more auxiliary TCL files into the Install Folder ***
REM if errorlevel 1 goto AuxTCLFailed
REM echo Installing auxiliary TCL Files: Completed OK!!

REM *** Copy the primary Netgen executable file into the Install Folder ***
echo Installing %PROJ_EXEC% into %INSTALL_FOLDER%\bin ....
if /i "%BUILD_ARCH%" == "win32" (
   xcopy "%PROJ_DIR%%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%\bin\" /i /d /y
   if errorlevel 1 goto ExecInstallFailed
)
if /i "%BUILD_ARCH%" == "x64" (
   xcopy "%PROJ_DIR%%BUILD_ARCH%\%BUILD_TYPE%\%PROJ_EXEC%" "%INSTALL_FOLDER%\bin\" /i /d /y
   if errorlevel 1 goto ExecInstallFailed
)   
echo Installing %PROJ_EXEC%: Completed OK!!

REM *** Copy the primary Netgen library and include files into the Install Folder ***
echo Installing %LIB_NAME%.lib into %INSTALL_FOLDER%\lib ....
if /i "%BUILD_ARCH%" == "win32" (
   xcopy "%PROJ_DIR%%BUILD_TYPE%\%LIB_NAME%.lib" "%INSTALL_FOLDER%\lib\" /i /d /y
   if errorlevel 1 goto LibInstallFailed
)
if /i "%BUILD_ARCH%" == "x64" (
   xcopy "%PROJ_DIR%%BUILD_ARCH%\%BUILD_TYPE%\%LIB_NAME%.lib" "%INSTALL_FOLDER%\lib\" /i /d /y
   if errorlevel 1 goto LibInstallFailed
)   
echo Installing %LIB_NAME%.lib: Completed OK!!

echo Installing %LIB_NAME%.h into %INSTALL_FOLDER%\include ....
xcopy "%NETGEN_LIBINC%\%LIB_NAME%.h" "%INSTALL_FOLDER%\include\" /i /d /y
if errorlevel 1 goto LibInstallFailed
echo Installing %LIB_NAME%.h: Completed OK!!

echo Installing NgSolve dependent header files into %INSTALL_FOLDER%\include ....
xcopy "%NETGEN_NGSINC%\include\nginterface_v2.hpp" "%INSTALL_FOLDER%\include\" /i /d /y
xcopy "%NETGEN_NGSINC%\general\dynamicmem.hpp" "%INSTALL_FOLDER%\include\" /i /d /y
xcopy "%NETGEN_NGSINC%\general\ngexception.hpp" "%INSTALL_FOLDER%\include\" /i /d /y
xcopy "%NETGEN_NGSINC%\visualization\soldata.hpp" "%INSTALL_FOLDER%\include\" /i /d /y

echo Installing external dependencies
if /i "%BUILD_ARCH%" == "x64" (   
   xcopy "%PROJ_DIR%..\..\ext_libs\pthreads-Win32\dll\x64\pthreadVC2.dll" "%INSTALL_FOLDER%\bin" /i /d /y   
   xcopy "%PROJ_DIR%..\..\ext_libs\zlib\x64\lib\zlib1.dll" "%INSTALL_FOLDER%\bin" /i /d /y   
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\bin\x64\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x64\tcl8.5" "%INSTALL_FOLDER%\lib\tcl8.5" /i /d /y /e /q
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x64\tix8.4.3" "%INSTALL_FOLDER%\lib\tix8.4.3" /i /d /y /e /q
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x64\tk8.5" "%INSTALL_FOLDER%\lib\tk8.5" /i /d /y /e /q
   copy "%PROJ_DIR%..\..\ext_libs\boost_1_56_0\lib\x64\boost_python-vc120-mt-1_56.dll" "%INSTALL_FOLDER%\bin\boost_python-vc120-mt-1_56.dll"
   REM if errorlevel 1 goto externalInstallFailed
)   

if /i "%BUILD_ARCH%" == "win32" (   
   xcopy "%PROJ_DIR%..\..\ext_libs\pthreads-Win32\dll\x86\pthreadVC2.dll" "%INSTALL_FOLDER%\bin" /i /d /y   
   xcopy "%PROJ_DIR%..\..\ext_libs\zlib\x86\lib\zlib1.dll" "%INSTALL_FOLDER%\bin" /i /d /y   
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\bin\x86\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x86\tcl8.5" "%INSTALL_FOLDER%\lib\tcl8.5" /i /d /y /e /q
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x86\tix8.4.3" "%INSTALL_FOLDER%\lib\tix8.4.3" /i /d /y /e /q
   xcopy "%PROJ_DIR%..\..\ext_libs\tcl\lib\x86\tk8.5" "%INSTALL_FOLDER%\lib\tk8.5" /i /d /y /e /q
   copy "%PROJ_DIR%..\..\ext_libs\boost_1_56_0\lib\Win32\boost_python-vc120-mt-1_56.dll" "%INSTALL_FOLDER%\bin\boost_python-vc120-mt-1_56.dll"
   REM if errorlevel 1 goto externalInstallFailed
)

REM *** Finally copy the FFMPEG dlls to the bin folder ***
if /i "%BUILD_ARCH%"=="x64" (
		echo Copying FFMPEG dlls
		C:\Windows\System32\xcopy "%PROJ_DIR%..\..\ext_libs\FFMPEG\dll\x64\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y		
)
if /i "%BUILD_ARCH%"=="win32" (
		echo Copying FFMPEG dlls
		C:\Windows\System32\xcopy "%PROJ_DIR%..\..\ext_libs\FFMPEG\dll\Win32\*.dll" "%INSTALL_FOLDER%\bin" /i /d /y		
)

echo Copying tutorials
xcopy "%PROJ_DIR%..\tutorials" "%INSTALL_FOLDER%\tutorials\" /i /d /y /e



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
:externalInstallFailed
echo POSTBUILD Script for %PROJ_NAME% FAILED..... Error copying external dependencies to install folder!!!
exit 1

:BuildEventOK
echo POSTBUILD Script for %PROJ_NAME% completed OK.....!!

