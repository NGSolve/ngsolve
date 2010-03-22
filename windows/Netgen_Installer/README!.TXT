*******************************************************************
* Netgen Mesher NSIS Installer
* Target Operating Systems: Windows XP / Vista / 7
* Target Achitectures: 32-bit / 64-bit
* NSIS Version: 2.46 (http://nsis.sourceforge.net/Main_Page)
*
* Author: Philippose Rajan
* Initial Release Date: 21 March 2010
*
*******************************************************************
* Description:
   This folder contains files necessary for generating and deploying 
the Windows Installer for the Netgen Mesher created using the 
Nullsoft Scriptable Install System (NSIS).

   The NSIS configuration/specification file of the installer 
supports automated generation of the NSIS Installer for both, the 
32-bit (Win32) as well as the 64-bit x64) variants of the Netgen 
Mesher.

   It should be noted however, that the Installer executable created 
by the NSIS compiler (for example, Netgen-4.9.12_Win32.exe) is itself, 
always a 32-bit executable. 

*******************************************************************
* Folder contents:
   ** netgen_installer.nsi     
         --> NSIS configuration/specification file used by the NSIS 
         compiler for generating the Installer executables.
         
   ** netgen_installer.dat
         --> Common Data File for the following batch routines. This 
         file contains user modifiable settings such as the version 
         number of Netgen to be used, the location of the sources and 
         binaries, etc...etc.         
         
   ** Netgen_Installer_All.bat
         --> Windows Command-line batch file for creating both, the 
         Installer for the 32-bit as well as the 64-bit variants of 
         the Netgen Mesher.
   
   ** Netgen_Installer_Win32.bat
         --> Windows Command-line batch file for creating the Installer 
         for the 32-bit (Win32) variant of the Netgen Mesher.
         
   ** Netgen_Installer_x64.bat
         --> Windows Command-line batch file for creating the Installer 
         for the 64-bit (x64) variant of the Netgen Mesher.
         
   ** LGPL_Licenses [Folder]
         --> This folder contains the LGPL License text in the various 
         languages supported by the Installer.
         --> Currently supported languages: 
             * English
             * German
             * French
             
*******************************************************************
* Netgen Version, Architecture and Folder structure Setup:
   The following variables need to be set in the automated compile 
batch files (Netgen_Installer_<All/Win32/x64>.bat) via the Common 
Data File (netgen_installer.dat):
   * VERSION - The current version of the Netgen Mesher (Ex. 4.9.12)
               for which the installers are being generated
   
   * ARCH    - The target architecture for which the installer is 
               to be generated
             - NOTE: This value is set within the BATCH Files, the 
               user must not make any modifications to this variable
             - Valid values: Win32 | x64
   
   * NETGEN  - Root of the Netgen installation on the computer on 
               which the NSIS compiler is being run
             - This is required to find the source and binary files 
               of the various variants of the Netgen Mesher during 
               generation of the Installer

   * NGSOLVE - Root of the NGSolve installation on the computer on 
               which the NSIS compiler is being run
             - This is required to find the source and binary files 
               of the various variants of the NGSolve add-on during 
               generation of the Installer
               
   * OCC     - Root of the OpenCascade installation on the computer 
               on which the NSIS compiler is being run
             - This is required to find the OpenCascade DLL files 
               which also need to be packahed into the Installer for 
               the Netgen Mesher in order to provide support for 
               OpenCascade Geometry

*******************************************************************   