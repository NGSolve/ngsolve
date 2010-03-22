@echo off
REM *************************************************************
REM Filename: Netgen_Installer_All.bat
REM
REM Automated NSIS Installer Compilation for creating the 
REM Graphical Windows Installer for the Netgen Meshing Software
REM
REM NOTE: This Batch file automatically generates Installers 
REM       for both, the 32-bit and the 64-bit versions of 
REM       Netgen
REM
REM Author: Philippose Rajan
REM Date: 11 March 2010
REM *************************************************************
SetLocal EnableDelayedExpansion

REM *** Name of the NSI File to be compiled by the NSIS compiler ***
set NSI_FILENAME=netgen_installer.nsi

REM ******* Read in the specification file "netgen_installer.dat" *******
REM * This file contains common settings for the 32-bit / 64-bit and 
REM * the combined versions of the automated compile batch files
REM * If the file does not exist, of if the required variables are not 
REM * found, default values are set
REM *********************************************************************
set DATA_FILE=netgen_installer.dat

if exist "%DATA_FILE%" (
   echo.
   echo Reading the Common Data File [%DATA_FILE%]....
   for /f "eol=# tokens=1,2 delims==" %%A in (netgen_installer.dat) do (
      if %%A==VERSION (
         set VERSION=%%B
         echo Found entry VERSION.... VERSION=!VERSION!
      )
      if %%A==NETGEN (
         set NETGEN=%%B
         echo Found entry NETGEN.... NETGEN=!NETGEN!         
      )
      if %%A==NGSOLVE (
         set NGSOLVE=%%B
         echo Found entry NGSOLVE.... NGSOLVE=!NGSOLVE!         
      )
      if %%A==OCC (
         set OCC=%%B
         echo Found entry OCC.... OCC=!OCC!         
      )      
   )
) else (
   echo.
   echo WARNING: Common Data File [%DATA_FILE%] not found.... Reverting to defaults!!
)   

if not defined VERSION set VERSION=4.9.XX
if not defined NETGEN set NETGEN=D:\netgenWin\05_Netgen_Main
if not defined NGSOLVE set NGSOLVE=D:\netgenWin\07_NGSolve_Main
if not defined OCC set OCC=D:\occ6.3.0
REM ********************************************************************   

echo.
echo 'Creating the NSIS Installer for Netgen 32-bit (Win32)....'
set ARCH=Win32
D:\NSIS\makensis.exe /DNETGEN_ARCH=%ARCH% /DNETGEN_VER=!VERSION! /DNETGEN_ROOT=!NETGEN! /DNGSOLVE_ROOT=!NGSOLVE! /DOCC_ROOT=!OCC! /Onetgen_%ARCH%_NSIS.log %NSI_FILENAME%
echo.
echo 'Done.... Please read the netgen_%ARCH%_NSIS.log file to check for errors....!'

echo.
echo 'Creating the NSIS Installer for Netgen 64-bit (Win64)....'
set ARCH=x64
D:\NSIS\makensis.exe /DNETGEN_ARCH=%ARCH% /DNETGEN_VER=!VERSION! /DNETGEN_ROOT=!NETGEN! /DNGSOLVE_ROOT=!NGSOLVE! /DOCC_ROOT=!OCC! /Onetgen_%ARCH%_NSIS.log %NSI_FILENAME%
echo.
echo 'Done.... Please read the netgen_%ARCH%_NSIS.log file to check for errors....!'
echo.
