; -----------------------------------------------
/*
   Nullsoft Installer for Netgen 32-bit and 64-bit
   
   Filename: netgen_installer.nsi
   Author  : Philippose Rajan
   Date    : 11 March 2010
   
   Description: This installer automates the process 
   of installing Netgen / Nglib and the associated 
   libraries on a Windows system
   ** Also gives the option of installing NgSolve 
   and samples of the Plug-In architecture
   ** NOTE: This installer can create both, the 32-bit 
   and the 64-bit versions of Netgen, depending on the 
   variable "NETGEN_ARCH"
      NETGEN_ARCH <Win32 | x64>

*/
; ------------------------------------------------



; ---------- Set the Compression Alg ----
SetCompressor /SOLID /FINAL lzma
; ---------------------------------------



; ---------- Include files --------------
; Include for sending Windows Messages
!include WinMessages.nsh

; Include the Modern UI subsystem
!include MUI2.nsh

; Include the Logic Library
!include LogicLib.nsh

; Miscellaneous string functions
;!include StrContains.nsh
; ---------------------------------------



; ---------- Variable declarations ------
Var foundTclTkTix
Var foundTogl
Var foundPThds
Var foundOCC
Var foundLapackBlas
Var foundMKLLibs
Var StartMenuFolder
Var NetgenUserDir
; ---------------------------------------



; ---------- House-keeping --------------
; Define the terms to be used for the various architectures 
; supported by Netgen
!define WIN32_ARCH   "Win32"
!define X64_ARCH     "x64"


; ########## MODIFY AS REQUIRED ###############
; Tell NSIS which architecture to use to create 
; the installer
; NOTE: Currently this variable is set at the 
; Command line when called "makensis.exe" using
; /DNETGEN_ARCH=<architecture>
;!define NETGEN_ARCH  ${WIN32_ARCH}
;!define NETGEN_ARCH  ${X64_ARCH}

; Change the version number of Netgen here
; NOTE: Currently this variable is set at the 
; Command line when called "makensis.exe" using
; /DNETGEN_VER=<version>
;!define NETGEN_VER   "4.9.12"



; Check if the Architecture and the version of Netgen have been specified
!ifndef NETGEN_ARCH
   !error "Error: NETGEN_ARCH not defined....Please define this on the command line.....Aborting!"
!endif
!ifndef NETGEN_VER
   !error "Error: NETGEN_VER not defined....Please define this on the command line.....Aborting!"
!endif


; Display the currently selected Netgen Architecture
!echo "Current Netgen Architecture: ${NETGEN_ARCH}"

; Display the current version of Netgen
!echo "Current Netgen Version: ${NETGEN_VER}"
; #############################################



; Part to be added to the "tcl" folders for 64-bit architecture
!if ${NETGEN_ARCH} == ${WIN32_ARCH}
   !define TCL_ARCH_EXT   ""
!else if ${NETGEN_ARCH} == ${X64_ARCH}
   !define TCL_ARCH_EXT   "-64"
!else
   !error "Error: Netgen Architecture not specified.....Aborting!"
!endif

; Part to be added to the "pthread" DLL for 32-bit or 64-bit architecture
!if ${NETGEN_ARCH} == ${WIN32_ARCH}
   !define PTHDS_ARCH_DIR "-w32"
   !define PTHDS_ARCH_EXT ""
!else if ${NETGEN_ARCH} == ${X64_ARCH}   
   !define PTHDS_ARCH_DIR "-w64"
   !define PTHDS_ARCH_EXT "_64"
!else
   !error "Error: Netgen Architecture not specified.....Aborting!"
!endif   

; Part to be added to the "Intel Math Kernel Libraries" DLL for 32-bit or 64-bit architecture
!if ${NETGEN_ARCH} == ${WIN32_ARCH}
   !define MKL_ARCH_DIR "-w32"
!else if ${NETGEN_ARCH} == ${X64_ARCH}   
   !define MKL_ARCH_DIR "-w64"
!else
   !error "Error: Netgen Architecture not specified.....Aborting!"
!endif

; Location of the OCC libraries for different architectures
!if ${NETGEN_ARCH} == ${WIN32_ARCH}
   !define OCC_ARCH "win32"
!else if ${NETGEN_ARCH} == ${X64_ARCH}   
   !define OCC_ARCH "win64"
!else
   !error "Error: Netgen Architecture not specified.....Aborting!"
!endif  



; ########## MODIFY AS REQUIRED ###############    
; Set up paths for access to various files
; Note #1: Change this according to the setup on 
; the machine where the installer is being compiled
; Note #2: These paths can be set from the command 
; line when calling "makensis.exe" by using:
; /DNETGEN_ROOT=<netgen root folder>
; /DNGSOLVE_ROOT=<ngsolve root folder>
; /DOCC_ROOT=<OpenCascade root folder>
!ifndef NETGEN_ROOT
   !define NETGEN_ROOT        "D:\netgenWin\05_Netgen_Main"
!endif
!ifndef NGSOLVE_ROOT
   !define NGSOLVE_ROOT       "D:\netgenWin\07_NGSolve_Main"
!endif
!ifndef OCC_ROOT
   !define OCC_ROOT           "D:\occ6.3.0"
!endif

; Display the currently selected Root Paths
!echo "Netgen Root: ${NETGEN_ROOT}"
!echo "NGSolve Root: ${NGSOLVE_ROOT}"
!echo "OpenCascade Root: ${OCC_ROOT}"



; Setup the directory structures that NSIS uses to access 
; the various flavours of Netgen, Nglib and NgSolve 
; Namely: 
; ** With and Without OpenCascade Support
; ** 32-bit or 64-bit versions
!define NETGEN_SRC         "${NETGEN_ROOT}\netgen"
!define NETGEN_BASIC       "${NETGEN_ROOT}\netgen-instNoOCC_${NETGEN_ARCH}"
!define NETGEN_OCC         "${NETGEN_ROOT}\netgen-instOCC_${NETGEN_ARCH}"

!define EXTLIBS_TCLTKTIX   "${NETGEN_ROOT}\ext_libs\tcl${TCL_ARCH_EXT}"
!define EXTLIBS_TOGL       "${NETGEN_ROOT}\ext_libs\tcl${TCL_ARCH_EXT}"
 
!define EXTLIBS_PTHDS      "${NETGEN_ROOT}\ext_libs\pthread${PTHDS_ARCH_DIR}"

!define NGLIB_SRC          "${NETGEN_ROOT}\netgen\nglib"
!define NGLIB_BIN_BASIC    "${NETGEN_ROOT}\nglib-instNoOCC_${NETGEN_ARCH}"
!define NGLIB_BIN_OCC      "${NETGEN_ROOT}\nglib-instOCC_${NETGEN_ARCH}"

!define DEMOAPP_SRC        "${NETGEN_ROOT}\demoapp"
!define DEMOAPP_BIN        "${NETGEN_ROOT}\netgen-demoapp_${NETGEN_ARCH}"

!define NGSOLVE_SRC        "${NGSOLVE_ROOT}\ngsolve"
!define NGSOLVE_BIN        "${NGSOLVE_ROOT}\ngsolve-inst_${NETGEN_ARCH}"
!define NGSOLVE_LAPACKLIBS "${NGSOLVE_ROOT}\ext_libs\lapack-blas"
!define NGSOLVE_MKLLIBS    "${NGSOLVE_ROOT}\ext_libs\mkl${MKL_ARCH_DIR}"

!define EXTLIBS_OCC        "${OCC_ROOT}\ros\${OCC_ARCH}\bin"

; Folder which contains the LGPL License in various languages
!define LGPL_LICENSE_SRC   ".\LGPL_Licenses"

; Location of the installer binary created by the NSIS compiler
!define OUTPUT_DIR         ".\Binaries"
; #############################################





; ########## USUALLY NOTHING NEEDS TO BE CHANGED BEYOND THIS POINT ###########

; Name of environment variable required by 
; Netgen for finding the required scripts
!define NETGENDIR       "NETGENDIR"
!define NETGEN_USER_DIR "NETGEN_USER_DIR"

; Paths within the Windows Registry which are 
; required for environment variables
!define HKLM_ENV           'HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"'
!define HKCU_ENV           'HKCU "Environment"'
!define SHELLFOLDERS       "Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders"

; Name to be used within the installer to refer to 
; the program
Name "Netgen Mesher (${NETGEN_ARCH})"

; Set the location of the output directory where the installer will be generated
!system 'md "${OUTPUT_DIR}"'

; Name of the generated installer
OutFile "${OUTPUT_DIR}\Netgen-${NETGEN_VER}_${NETGEN_ARCH}.exe"

; Set the style 
XPStyle on

; Default installation folder
InstallDir "$PROGRAMFILES\Netgen-${NETGEN_VER}_${NETGEN_ARCH}"

; Try to get the current installed folder (if any)
InstallDirRegKey HKCU "Software\Netgen\${NETGEN_VER}_${NETGEN_ARCH}" "Install_Dir"

; Request application privileges for Windows Vista
RequestExecutionLevel user
; ---------------------------------------



; ---------- Interface settings ---------
; Currently uses the Modern UI v2 User Interface

; Define the bitmap to use for the welcome and finish pages
!define MUI_WELCOMEFINISHPAGE_BITMAP   "${NSISDIR}\Contrib\Graphics\Wizard\orange.bmp"
!define MUI_UNWELCOMEFINISHPAGE_BITMAP "${NSISDIR}\Contrib\Graphics\Wizard\orange.bmp"

; Display a header bitmap image on top right of each page
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_RIGHT
!define MUI_HEADERIMAGE_BITMAP   "${NSISDIR}\Contrib\Graphics\Header\orange-r.bmp"
!define MUI_HEADERIMAGE_UNBITMAP "${NSISDIR}\Contrib\Graphics\Header\orange-r.bmp"

!define MUI_ABORTWARNING

; Select the Icon to use for the generated installer
!define MUI_ICON   "${NSISDIR}\Contrib\Graphics\Icons\orange-install-nsis.ico"
!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\orange-uninstall-nsis.ico"
; ---------------------------------------


 
; ---------- Pages to display -----------
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "$(license)"
!insertmacro MUI_PAGE_COMPONENTS

; Select the Installation Directory
!define MUI_DIRECTORYPAGE_VARIABLE $INSTDIR
!define MUI_DIRECTORYPAGE_TEXT_DESTINATION $(STR_InstallFolder)
!insertmacro MUI_PAGE_DIRECTORY

; Select the Netgen user Directory
!define MUI_PAGE_HEADER_TEXT $(STR_UserPageHeaderText)
!define MUI_PAGE_HEADER_SUBTEXT $(STR_UserPageHeaderSubText)
!define MUI_DIRECTORYPAGE_TEXT_TOP $(STR_UserPageText)
!define MUI_DIRECTORYPAGE_TEXT_DESTINATION $(STR_UserFolder)
!define MUI_DIRECTORYPAGE_VARIABLE $NetgenUserDir
!insertmacro MUI_PAGE_DIRECTORY

; Set up the Start Menu Page 
!define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU"
!define MUI_STARTMENUPAGE_REGISTRY_KEY  "Software\Netgen\${NETGEN_VER}_${NETGEN_ARCH}"
!define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start_Menu_Folder"
!define MUI_STARTMENUPAGE_DEFAULTFOLDER "Netgen"
!insertmacro MUI_PAGE_STARTMENU NetgenStartMenu $StartMenuFolder

!insertmacro MUI_PAGE_INSTFILES

!define MUI_FINISHPAGE_RUN "$INSTDIR\bin\netgen.exe"
!define MUI_FINISHPAGE_RUN_TEXT $(STR_Finish_RunText)
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH
; ---------------------------------------



; ---------- Languages to use -----------
!insertmacro MUI_LANGUAGE "English"
!insertmacro MUI_LANGUAGE "German"
!insertmacro MUI_LANGUAGE "French"
; ---------------------------------------



; ---------- Reserve Files --------------
!insertmacro MUI_RESERVEFILE_LANGDLL
; ---------------------------------------


; ---------- License files --------------
LicenseLangString license ${LANG_ENGLISH} ${LGPL_LICENSE_SRC}\LGPL_English.txt
LicenseLangString license ${LANG_GERMAN} ${LGPL_LICENSE_SRC}\LGPL_German.txt
LicenseLangString license ${LANG_FRENCH} ${LGPL_LICENSE_SRC}\LGPL_French.txt
LicenseData $(license)
; ---------------------------------------



; ---------- File Copy Macros -----------
; Macro to copy basic Netgen files into the 
; install location
; These files are common to the versions of 
; netgen with and without OCC support
!macro copyNetgenBasics reqDir
    ; Set the output path to the install 
	; folder path
    SetOutPath $INSTDIR
	
	; Create the directory structure
	; used for the installation
	CreateDirectory "$INSTDIR\bin"
	CreateDirectory "$INSTDIR\include"
	CreateDirectory "$INSTDIR\lib"
	CreateDirectory "$INSTDIR\doc"

	; Create the user specific directory where netgen 
	; stores configuration information	
	CreateDirectory "$NetgenUserDir"
	  
    ; Copy the files into the right 
	; folders
	SetOutPath "$INSTDIR\bin"
    File ${reqDir}\bin\*.tcl
	File ${reqDir}\bin\*.exe
	File ${reqDir}\bin\*.ocf
	  
	SetOutPath "$INSTDIR\include"
	File ${reqDir}\include\*.h
	File ${reqDir}\include\*.hpp
	  
	SetOutPath "$INSTDIR\lib"
	File ${reqDir}\lib\*.lib
	  
	SetOutPath "$INSTDIR\doc"
	File ${NETGEN_SRC}\doc\*.pdf   
	
	SetOutPath "$INSTDIR\tutorials"
	File ${NETGEN_SRC}\tutorials\*.*
!macroend
; ---------------------------------------



; Macro to copy the common auxiliary  
; libraries
; The main libraries handled here are the 
; Tcl, Tk, Tix, Togl and pthread-w32 libraries
!macro copyAuxLibBasics
	; Install the required Auxiliary Libraries
	; (If not already available and accessible)
	${If} $foundTclTkTix != "1"
	    SetOutPath "$INSTDIR\bin"
	    File ${EXTLIBS_TCLTKTIX}\bin\tcl85.dll
	    File ${EXTLIBS_TCLTKTIX}\bin\tclpip85.dll
	    File ${EXTLIBS_TCLTKTIX}\bin\tk85.dll
	    File ${EXTLIBS_TCLTKTIX}\bin\tix84.dll
		
		; Copy the basic Tcl headers and lib files to help
		; with compilation of NgSolve
		SetOutPath "$INSTDIR\include"
		File ${EXTLIBS_TCLTKTIX}\include\tcl.h
		File ${EXTLIBS_TCLTKTIX}\include\tclDecls.h
		File ${EXTLIBS_TCLTKTIX}\include\tclPlatDecls.h
		
	    SetOutPath "$INSTDIR\lib"
		File ${EXTLIBS_TCLTKTIX}\lib\tcl85.lib		
	    File /r /x *.lib /x demos /x tzdata ${EXTLIBS_TCLTKTIX}\lib\*.*
	${EndIf}	

	${If} $foundTogl != "1"
        SetOutPath "$INSTDIR\bin"
	    File ${EXTLIBS_TOGL}\bin\Togl17.dll
    ${EndIf}		
	
    ${If} $foundPThds != "1"
	    SetOutPath "$INSTDIR\bin"
	        File ${EXTLIBS_PTHDS}\lib\pthreadVC2${PTHDS_ARCH_EXT}.dll
    ${EndIf}		
!macroend
; ---------------------------------------



; Macro to write the install folder and the 
; NETGENDIR environment variable into the registry, 
; and write the uninstaller.exe file
!macro registerInstall
    ; Set the NETGENDIR environment variable to the current 
    ; installation bin folder
    WriteRegExpandStr ${HKCU_ENV} "${NETGENDIR}" "$INSTDIR\bin"
	
	; Set the NETGEN_USER_DIR environment variable to the user 
	; selected folder
	WriteRegExpandStr ${HKCU_ENV} "${NETGEN_USER_DIR}" $NetgenUserDir

	; Broadcast the changes in the environment variables to the rest of the system
	SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
	
    ; Write the installation path into the registry
    WriteRegStr HKCU "Software\Netgen\${NETGEN_VER}_${NETGEN_ARCH}" "Install_Dir" "$INSTDIR"
    WriteRegStr HKCU "Software\Netgen\${NETGEN_VER}_${NETGEN_ARCH}" "User_Dir" "$NetgenUserDir"	
	
	; Write the un-installer
	WriteUninstaller "$INSTDIR\uninstall.exe"
	
	SetOutPath "$INSTDIR\tutorials"
	
	; Write the Start Menu items
	!insertmacro MUI_STARTMENU_WRITE_BEGIN NetgenStartMenu
	  CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
	  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Netgen ${NETGEN_VER}.lnk" "$INSTDIR\bin\netgen.exe"
	  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Netgen User Manual.lnk" "$INSTDIR\doc\ng4.pdf"
	  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk" "$INSTDIR\uninstall.exe"
	!insertmacro MUI_STARTMENU_WRITE_END  	
!macroend
; ---------------------------------------
; ---------------------------------------



; ---------- Main Section ---------------
SectionGroup /e "Netgen Executable"

    Section "Netgen Basic" secNetgen
   
        LogSet on
   
        ; Copy the basic Netgen files into the 
        ; install location	  
        !insertmacro copyNetgenBasics ${NETGEN_BASIC}

	    ; Check if required libraries are present
	    ; if not, mark them to be installed too
	    Call checkTclTkTix
	    Call checkTogl
	    Call checkPThds
		
	    ; Copy the basic auxiliary libraries 
	    ; (if required) to the install location
        !insertmacro copyAuxLibBasics
	
		; register the install by adding registry 
		; keys, and write the uninstaller
		!insertmacro registerInstall
    SectionEnd

	Section /o "Netgen with OCC" secNetgenOCC

        LogSet on	
	
		; Copy the basic Netgen files into the 
		; install location	 
		!insertmacro copyNetgenBasics ${NETGEN_OCC}

		; Check if required libraries are present
		; if not, mark them to be installed too
		Call checkTclTkTix
		Call checkTogl
		Call checkPThds
		Call checkOCC
		
		; Copy the basic auxiliary libraries 
		; (if required) to the install location
		!insertmacro copyAuxLibBasics

		StrCmp $foundOCC "1" OCCInstalled 0
			SetOutPath "$INSTDIR\bin"
			File ${EXTLIBS_OCC}\*.dll		 
	OCCInstalled: 	  

		; register the install by adding registry 
		; keys, and write the uninstaller
		!insertmacro registerInstall
	SectionEnd
	
SectionGroupEnd



SectionGroup "Netgen API DLL (Nglib)"   
	Section /o "Nglib" secNglib
	
	    ; Copy the nglib Basic DLL to the bin folder
		; in the install location
		SetOutPath "$INSTDIR\bin"
		File ${NGLIB_BIN_BASIC}\bin\nglib.dll
		
		; Copy the nglib LIB file to the lib folder 
		; in the install location
		SetOutPath "$INSTDIR\lib"
		File ${NGLIB_BIN_BASIC}\lib\nglib.lib

		; Copy the nglib include header file to the 
		; include folder in the install location
		SetOutPath "$INSTDIR\include"
		File ${NGLIB_SRC}\nglib.h
	SectionEnd
    
    Section /o "Nglib with OCC" secNglibOCC
    
	    ; Copy the nglib DLL version with OCC support to the bin folder
		; in the install location
		SetOutPath "$INSTDIR\bin"
		File ${NGLIB_BIN_OCC}\bin\nglib.dll
		
		; Copy the nglib LIB file to the lib folder 
		; in the install location
		SetOutPath "$INSTDIR\lib"
		File ${NGLIB_BIN_OCC}\lib\nglib.lib

		; Copy the nglib include header file to the 
		; include folder in the install location
		SetOutPath "$INSTDIR\include"
		File ${NGLIB_SRC}\nglib.h        
    SectionEnd
	
	Section /o "Nglib examples" secNglibEx
	
	SectionEnd
	
SectionGroupEnd



SectionGroup "Add-On Applications"
	
	Section /o "NgSolve" secNgSolve
	
	    ; Check if required libraries are present
	    ; if not, mark them to be installed too		
		Call checkLapackBlas
        Call checkMKLLibs
	
		; Copy the DLL and Tcl files to the 
		; Netgen bin folder
		SetOutPath "$INSTDIR\bin"
		File ${NGSOLVE_BIN}\bin\*.dll
		File ${NGSOLVE_BIN}\bin\*.tcl
		
		${If} ${NETGEN_ARCH} != ${X64_ARCH}
            ${If} $foundLapackBlas != "1"
	            SetOutPath "$INSTDIR\bin"
	            File ${NGSOLVE_LAPACKLIBS}\bin\*.dll		 
            ${EndIf}
		${EndIf}
        
        ${If} $foundMKLLibs != "1"
           SetOutPath "$INSTDIR\bin"
           File ${NGSOLVE_MKLLIBS}\bin\*.dll
        ${EndIf}   
		
        ; Copy the NGSolve LIB file to the 
        ; Netgen lib folder
        SetOutPath "$INSTDIR\lib"
        File ${NGSOLVE_BIN}\lib\*.lib
        
        ; Copy the NGSolve header file structure  
        ; to the Netgen include folder
        SetOutPath "$INSTDIR\include\include"
        File ${NGSOLVE_SRC}\include\*.hpp
        
        SetOutPath "$INSTDIR\include\solve"
        File ${NGSOLVE_SRC}\solve\*.hpp
        
        SetOutPath "$INSTDIR\include\parallel"
        File ${NGSOLVE_SRC}\parallel\*.hpp
        
        SetOutPath "$INSTDIR\include\ngstd"
        File ${NGSOLVE_SRC}\ngstd\*.hpp
        
        SetOutPath "$INSTDIR\include\multigrid"
        File ${NGSOLVE_SRC}\multigrid\*.hpp
        
        SetOutPath "$INSTDIR\include\linalg"
        File ${NGSOLVE_SRC}\linalg\*.hpp
        
        SetOutPath "$INSTDIR\include\fem"
        File ${NGSOLVE_SRC}\fem\*.hpp
        
        SetOutPath "$INSTDIR\include\comp"
        File ${NGSOLVE_SRC}\comp\*.hpp
        
        SetOutPath "$INSTDIR\include\basiclinalg"
        File ${NGSOLVE_SRC}\basiclinalg\*.hpp
        
       
	    ; Copy the NgSolve pde tutorials to the 
		; installation folder
		SetOutPath "$INSTDIR\pde_tutorials"
		File ${NGSOLVE_SRC}\pde_tutorial\*.*
	SectionEnd
	
	
	Section /o "demoapp" secDemoApp
	
		; Copy the demoapp DLL and Tcl files 
		; to the Netgen bin folder
		SetOutPath "$INSTDIR\bin"
		File ${DEMOAPP_BIN}\bin\*.dll
		File ${DEMOAPP_BIN}\bin\*.tcl
		
		; Copy the source files to the plugins folder
		; under a subdirectory named "demoapp"
		SetOutPath "$INSTDIR\plugins\demoapp"
		File ${DEMOAPP_SRC}\*.cpp
		File ${DEMOAPP_SRC}\*.h
		File ${DEMOAPP_SRC}\*.tcl
		File ${DEMOAPP_SRC}\README
		
		SetOutPath "$INSTDIR\plugins\demoapp\windows"
		File ${DEMOAPP_SRC}\windows\*.sln
		File ${DEMOAPP_SRC}\windows\*.vcproj
		File ${DEMOAPP_SRC}\windows\*.bat
	SectionEnd
	
SectionGroupEnd   
; ---------------------------------------



; ---------- Un-Installer Section -------
Section Uninstall

	; The first step is to always delete the 
	; uninstaller
	Delete "$INSTDIR\uninstall.exe"
	
	; Remove the folders one by one...
	; Could be done in one shot, but prefer 
	; to be a little clearer here
	RMDir /r /REBOOTOK "$INSTDIR\bin"
	RMDir /r /REBOOTOK "$INSTDIR\lib"
	RMDir /r /REBOOTOK "$INSTDIR\include"
	RMDir /r /REBOOTOK "$INSTDIR\doc"
	RMDir /r /REBOOTOK "$INSTDIR\plugins"
	RMDir /r /REBOOTOK "$INSTDIR\tutorials"
	
	; Finally remove the install folder
	RMDir /r /REBOOTOK $INSTDIR
	
	; Remove the NETGENDIR registry key
	DeleteRegValue ${HKCU_ENV} "${NETGENDIR}"
	DeleteRegValue ${HKCU_ENV} "${NETGEN_USER_DIR}"
	
	!insertmacro MUI_STARTMENU_GETFOLDER NetgenStartMenu $StartMenuFolder
	RMDir /r /REBOOTOK "$SMPROGRAMS\$StartMenuFolder"
	
	; Remove the Netgen entry from the SOFTWARE registry
	DeleteRegKey HKCU "SOFTWARE\Netgen\${NETGEN_VER}_${NETGEN_ARCH}"
	DeleteRegKey /ifempty HKCU "SOFTWARE\Netgen"
    
	; Broadcast the changes in the environment variables to the rest of the system
	SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000    
SectionEnd	


; ---------- Section Descriptions -------

LangString DESC_secNetgen     ${LANG_ENGLISH} "Basic version of the Netgen Mesher"
LangString DESC_secNetgen     ${LANG_GERMAN}  "Basis Version des Netgen Vernetzungstools"
LangString DESC_secNetgen     ${LANG_FRENCH}  "Version de base de l'outil de maillage Netgen"

LangString DESC_secNetgenOCC  ${LANG_ENGLISH} "Netgen Mesher version with OpenCascade geometry support (STEP/IGES/BREP files)"
LangString DESC_secNetgenOCC  ${LANG_GERMAN}  "Der Netgen Vernetzungstool mit eingebauter Unterstützung für OpenCascade Geometrie (STEP/IGES/BREP Dateien)"
LangString DESC_secNetgenOCC  ${LANG_FRENCH}  "Netgen outil de maillage avec le soutien de la géométrie OpenCascade (STEP/IGES/BREP fichiers)"

LangString DESC_secNglib      ${LANG_ENGLISH} "Library (DLL) version of the Netgen Mesher containing the API for use in Third-party Software"
LangString DESC_secNglib      ${LANG_GERMAN}  "Bibliotheksversion (DLL) des Netgen Vernetzungstools für Anwendung von innerhalb externen Software"
LangString DESC_secNglib      ${LANG_FRENCH}  "Bibliothèque (DLL) version de l'outil de maillage Netgen contenant l'API pour une utilisation dans Logiciels tiers"

LangString DESC_secNglibOCC   ${LANG_ENGLISH} "Library (DLL) version of the Netgen Mesher with OpenCascade geometry support containing the API for use in Third-party Software"
LangString DESC_secNglibOCC   ${LANG_GERMAN}  "Bibliotheksversion (DLL) des Netgen Vernetzungstools mit eingebauter Unterstützung für OpenCascade Geometrie für Anwendung von innerhalb externen Software"
LangString DESC_secNglibOCC   ${LANG_FRENCH}  "Bibliothèque (DLL) version de l'outil de maillage Netgen contenant l'API pour une utilisation dans Logiciels tiers avec le soutien de la géométrie OpenCascade"


!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${secNetgen}    $(DESC_secNetgen)

    !insertmacro MUI_DESCRIPTION_TEXT ${secNetgenOCC} $(DESC_secNetgenOCC)

    !insertmacro MUI_DESCRIPTION_TEXT ${secNglib}     $(DESC_secNglib)
   
    !insertmacro MUI_DESCRIPTION_TEXT ${secNglibOCC}  $(DESC_secNglibOCC)   
	
	!insertmacro MUI_DESCRIPTION_TEXT ${secNgSolve} "Pre-compiled (DLL) version of the NgSolve module - The Netgen Solver Module"
	
	!insertmacro MUI_DESCRIPTION_TEXT ${secDemoApp} "Sample of a plugin application for customising and extending Netgen - Source files and pre-compiled DLL"
	
	!insertmacro MUI_DESCRIPTION_TEXT ${secNglibEx} "Sample applications for demonstrating the NgLib API DLL"
!insertmacro MUI_FUNCTION_DESCRIPTION_END


; ---------- Miscellaneous Description Strings -------

LangString STR_InstallFolder       ${LANG_ENGLISH} "Installation Directory"
LangString STR_InstallFolder       ${LANG_GERMAN}  "Installationsordner"
LangString STR_InstallFolder       ${LANG_FRENCH}  "Dossier d'installation"

LangString STR_UserFolder          ${LANG_ENGLISH} "Directory for User specific settings"
LangString STR_UserFolder          ${LANG_GERMAN}  "Ordner für Benutzerspezifische Einstellungen"
LangString STR_UserFolder          ${LANG_FRENCH}  "Dossier pour les paramètres spécifiques de l'utilisateur"

LangString STR_UserPageText        ${LANG_ENGLISH} "Netgen uses the directory specified here in order to save user-specific configuration information. Make sure that you have write access for the selected folder."
LangString STR_UserPageText        ${LANG_GERMAN}  "Netgen verwendet das hier angegebene Verzeichnis um anwenderspezifische Konfigurationsdatein zu speichern. Vergewissern Sie sich dass Sie Schreibrechte in diesem Verzeichnis haben"
LangString STR_UserPageText        ${LANG_FRENCH}  "Netgen utilise le dossier spécifié ici afin de sauver l'utilisateur des informations de configuration spécifiques. Assurez-vous que vous avez accès en écriture pour le dossier sélectionné."

LangString STR_UserPageHeaderText  ${LANG_ENGLISH} "Select Location for User Specific Settings"
LangString STR_UserPageHeaderText  ${LANG_GERMAN}  "Zielverzeichnis für Benutzerspezifische Einstellungen"
LangString STR_UserPageHeaderText  ${LANG_FRENCH}  "Sélectionnez le dossier contenant les paramètres spécifiques à l'utilisateur"

LangString STR_UserPageHeaderSubText  ${LANG_ENGLISH} "Select the folder which Netgen Mesher can use to save user specific data"
LangString STR_UserPageHeaderSubText  ${LANG_GERMAN}  "Wählen Sie das Verzeichnis aus, in das Netgen Mesher die Benutzerspezifische Daten speichern soll"
LangString STR_UserPageHeaderSubText  ${LANG_FRENCH}  "Sélectionnez le dossier dans lequel l'outil de maillage Netgen pouvez utiliser pour enregistrer des données utilisateur spécifiques"

LangString STR_Finish_RunText         ${LANG_ENGLISH} "Run Netgen-${NETGEN_VER} now"
LangString STR_Finish_RunText         ${LANG_GERMAN}  "Netgen-${NETGEN_VER} jetzt starten"
LangString STR_Finish_RunText         ${LANG_FRENCH}  "Exécutez le Netgen-${NETGEN_VER} maintenant"
; ---------------------------------------



; ----------- Language selection --------
Function .onInit
    ; Display the language select window
    !insertmacro MUI_LANGDLL_DISPLAY
   
    ; The base Netgen package (without OCC) is selected by default
    StrCpy $1 ${secNetgen}
    !insertmacro SetSectionFlag ${secNglibOCC} ${SF_RO}
    !insertmacro ClearSectionFlag ${secNglib} ${SF_RO}
	
	; Extract the user application data folder 
	ReadRegStr $3 HKCU "${SHELLFOLDERS}" AppData
	StrCmp $3 "" 0 +3
	  StrCpy $NetgenUserDir "$INSTDIR"
	  Goto End
	StrCpy $NetgenUserDir "$3\Netgen-${NETGEN_VER}"
End:	  
FunctionEnd

; The init function for the uninstaller
Function un.oninit
    ; Display the language select window
    !insertmacro MUI_LANGDLL_DISPLAY
FunctionEnd
; ---------------------------------------



; ---------- Package Selection ----------
Function .onSelChange
    !insertmacro StartRadioButtons $1
        !insertmacro RadioButton ${secNetgen}
        !insertmacro RadioButton ${secNetgenOCC}
    !insertmacro EndRadioButtons


    ; Select the right version of Nglib based on the current 
    ; Netgen package version selection
    !insertmacro SectionFlagIsSet ${secNetgen} ${SF_SELECTED} "NglibBasicSel" ""
    !insertmacro SectionFlagIsSet ${secNetgenOCC} ${SF_SELECTED} "NglibOCCSel" ""
        Goto NglibSelEnd
   NglibBasicSel:     
        !insertmacro UnselectSection ${secNglibOCC}
        !insertmacro SetSectionFlag ${secNglibOCC} ${SF_RO}
        !insertmacro ClearSectionFlag ${secNglib} ${SF_RO}
        Goto NglibSelEnd   
   NglibOCCSel:
        !insertmacro UnselectSection ${secNglib}
        !insertmacro SetSectionFlag ${secNglib} ${SF_RO}
        !insertmacro ClearSectionFlag ${secNglibOCC} ${SF_RO}        
   NglibSelEnd:       
    
FunctionEnd
; ---------------------------------------



; ------------- StrContains -------------
; StrContains
; This function does a case sensitive searches for an occurrence of a substring in a string. 
; It returns the substring if it is found. 
; Otherwise it returns null(""). 
; Written by kenglish_hi
; Adapted from StrReplace written by dandaman32
 
 
Var STR_HAYSTACK
Var STR_NEEDLE
Var STR_CONTAINS_VAR_1
Var STR_CONTAINS_VAR_2
Var STR_CONTAINS_VAR_3
Var STR_CONTAINS_VAR_4
Var STR_RETURN_VAR
 
Function StrContains
  Exch $STR_NEEDLE
  Exch 1
  Exch $STR_HAYSTACK
  ; Uncomment to debug
  ;MessageBox MB_OK 'STR_NEEDLE = $STR_NEEDLE STR_HAYSTACK = $STR_HAYSTACK '
    StrCpy $STR_RETURN_VAR ""
    StrCpy $STR_CONTAINS_VAR_1 -1
    StrLen $STR_CONTAINS_VAR_2 $STR_NEEDLE
    StrLen $STR_CONTAINS_VAR_4 $STR_HAYSTACK
    loop:
      IntOp $STR_CONTAINS_VAR_1 $STR_CONTAINS_VAR_1 + 1
      StrCpy $STR_CONTAINS_VAR_3 $STR_HAYSTACK $STR_CONTAINS_VAR_2 $STR_CONTAINS_VAR_1
      StrCmp $STR_CONTAINS_VAR_3 $STR_NEEDLE found
      StrCmp $STR_CONTAINS_VAR_1 $STR_CONTAINS_VAR_4 done
      Goto loop
    found:
      StrCpy $STR_RETURN_VAR $STR_NEEDLE
      Goto done
    done:
   Pop $STR_NEEDLE ;Prevent "invalid opcode" errors and keep the
   Exch $STR_RETURN_VAR  
FunctionEnd
 
!macro _StrContainsConstructor OUT NEEDLE HAYSTACK
  Push "${HAYSTACK}"
  Push "${NEEDLE}"
  Call StrContains
  Pop "${OUT}"
!macroend
 
!define StrContains '!insertmacro "_StrContainsConstructor"'
; ---------------------------------------




; ---------- Check for TCL/Tk DLLs ------
Function checkTclTkTix
	ClearErrors
	; Check if the Tix, Tk and Tcl DLLs are available in the PATH
	; Note: Modify version numbers if required
	SearchPath $2 tix84.dll
	SearchPath $2 tk85.dll
	SearchPath $2 tcl85.dll
	IfErrors TclTkTixErrors
	Push $2
	Call GetPathOnly
	Pop "$R0"
	; Check if the required auxiliary Tcl files are present 
	; in the location expected by Tcl/Tk/Tix
	IfFileExists $R0\..\lib\tcl8.5\init.tcl 0 TclTkTixErrors	
	IfFileExists $R0\..\lib\tk8.5\*.tcl 0 TclTkTixErrors
	IfFileExists $R0\..\lib\tix8.4.* 0 TclTkTixErrors	
	${StrContains} $2 "${TCL_ARCH_EXT}" $R0
	${If} $2 == "" 
	${AndIf} ${NETGEN_ARCH} == ${X64_ARCH}
	    Goto TclTkTixErrors
	${EndIf}
	${If} $2 != "" 
	${AndIf} ${NETGEN_ARCH} == ${WIN32_ARCH}
	    Goto TclTkTixErrors
	${EndIf}	
    StrCpy $foundTclTkTix "1"
	LogText "Tcl/Tk and Tix Packages already installed on system....skipping local installation"
	ClearErrors
	Return
TclTkTixErrors:	
	StrCpy $foundTclTkTix "0"
	LogText "Tcl/Tk and Tix Packages not found on system....installing locally"	
	ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Check for Togl DLLs --------
Function checkTogl
	ClearErrors
	SearchPath $2 Togl17.dll
	IfErrors ToglErrors
	Push $2
	Call GetPathOnly
	Pop "$R0"
	${StrContains} $2 "${TCL_ARCH_EXT}" $R0
	${If} $2 == "" 
	${AndIf} ${NETGEN_ARCH} == ${X64_ARCH}
	    Goto ToglErrors
	${EndIf}
	${If} $2 != "" 
	${AndIf} ${NETGEN_ARCH} == ${WIN32_ARCH}
	    Goto ToglErrors
	${EndIf}	
	StrCpy $foundTogl "1"
	LogText "Togl-1.7 already installed on system....skipping local installation"	
	ClearErrors
	Return
ToglErrors:		
	StrCpy $foundTogl "0"
	LogText "Togl-1.7 not found on system....installing locally"		
	ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Check for pthreads DLLs ----
Function checkPThds
	ClearErrors
	${If} ${NETGEN_ARCH} == ${X64_ARCH}
		SearchPath $2 pthreadVC2_64.dll
	${Else}	
	    SearchPath $2 pthreadVC2.dll
	${EndIf}
	IfErrors PThdsErrors
		StrCpy $foundPThds "1"
	    LogText "pthread${PTHDS_ARCH_DIR} already installed on system....skipping local installation"		
		ClearErrors
		Return
PThdsErrors:		
		StrCpy $foundPThds "0"
		LogText "pthread${PTHDS_ARCH_DIR} not found on system....installing locally"	
		ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Check for OCC DLLs ---------
; Note: Not all required OCC DLLs will be 
; checked for... only a few representative 
; ones.
Function checkOCC
	ClearErrors
	SearchPath $2 TKernel.dll
	SearchPath $2 TKLCAF.dll
	SearchPath $2 TKSTEP.dll
	SearchPath $2 TKIGES.dll
	SearchPath $2 TKBRep.dll
	IfErrors OCCErrors
	Push $2
	Call GetPathOnly
	Pop "$R0"	
	${StrContains} $2 ${OCC_ARCH} $R0
	${If} $2 == "" 
	    Goto OCCErrors
	${EndIf}	
		StrCpy $foundOCC "1"
		LogText "OpenCascade Libraries already installed on system....skipping local installation"	
		ClearErrors
		Return
OCCErrors:		
		StrCpy $foundOCC "0"
		LogText "OpenCascade Libraries not found on system....installing locally"			
		ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Check for Lapack/BLAS DLLs -
; Checks if the Lapack/BLAS DLLs are available 
; via the PATH environment variable
Function checkLapackBlas
	ClearErrors
	SearchPath $2 lapack_win32.dll
	SearchPath $2 blas_win32.dll
	IfErrors LapackBlasErrors
		StrCpy $foundLapackBlas "1"
		LogText "Lapack / BLAS Libraries already installed on system....skipping local installation"		
		ClearErrors
		Return
LapackBlasErrors:
		StrCpy $foundLapackBlas "0"
		LogText "Lapack / BLAS Libraries not found on system....installing locally"		
		ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Check for Intel Match Kernel Libraries (MKL) DLLs -
; Checks if the MKL DLLs are available 
; via the PATH environment variable
Function checkMKLLibs
	ClearErrors
	SearchPath $2 libiomp5md.dll
	SearchPath $2 mkl_core.dll
    SearchPath $2 mkl_def.dll
    SearchPath $2 mkl_intel_thread.dll
    SearchPath $2 msvcr71.dll
	IfErrors MKLErrors
		StrCpy $foundMKLLibs "1"
		LogText "Math Kernel Libraries (MKL) already installed on system....skipping local installation"		
		ClearErrors
		Return
MKLErrors:
		StrCpy $foundMKLLibs "0"
		LogText "Math Kernel Libraries (MKL) not found on system....installing locally"		
		ClearErrors
FunctionEnd
; ---------------------------------------



; ---------- Func1 - Get Path Only ------
Function GetPathOnly
  Exch $R0
  Push $R1
  Push $R2
  StrLen $R1 $R0
  IntOp $R1 $R1 + 1
  loop:
    IntOp $R1 $R1 - 1
    StrCpy $R2 $R0 1 $R1
    IntCmp $R1 1 exit2
    StrCmp $R2 "\" exit1 ; Look for the first backslash from the end of string
  Goto loop
  exit1:
    StrCpy $R0 $R0 $R1
  exit2:
    Pop $R2
    Pop $R1
    Exch $R0
FunctionEnd   
; ---------------------------------------

