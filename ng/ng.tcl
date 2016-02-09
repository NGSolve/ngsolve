lappend auto_path $env(NETGENDIR) 

set batchmode [Ng_GetCommandLineParameter batchmode]
if {$batchmode=="undefined"} {
    if {[catch {package require Tix } result ]} {
        puts "cannot load package Tix"
        puts "error : $result"
    }
    if {[catch {package require tkdnd } result ]} {
#        puts "cannot load package tkdnd"
#        puts "error : $result"
    }	
}
# if {[catch {package require Togl 2.0 } result ]} {
#    puts "cannot load package Togl 2.0"
#    puts "error : $result"
# }

# load [file dirname [info script]]/gears[info sharedlibextension]
# puts "load togl lib"
# load /Users/joachim/tcl_native3/Togl2.1/libTogl2.1.dylib
# puts "have togl lib"
# if {[catch {package require Togl 2.0 } result ]} {
#    puts "cannot load package Togl 2.0"
#    puts "error : $result"
# }



# userlevel 1..standard user 2..power-user 3..developer

set userlevel 3
if { [Ng_GetCommandLineParameter expert]=="defined" } {
    set userlevel 3
}

set progname "NETGEN"

set ngdir ""
if { [lsearch [array names env] NETGENDIR] != -1 } {
    set ngdir $env(NETGENDIR) 
}
if { [string length $ngdir] == 0 } {
    set ngdir "." 
}

set nguserdir ""
if { [lsearch [array names env] NETGEN_USER_DIR] != -1 } {
    set nguserdir $env(NETGEN_USER_DIR) 
}
if { [string length $nguserdir] == 0 } {
    set nguserdir "." 
}




set batchmode [Ng_GetCommandLineParameter batchmode]

set solvemode 0
if { [Ng_GetCommandLineParameter solve] != "undefined" || \
	 [Ng_GetCommandLineParameter recent] == "defined" } {
    set solvemode defined
}

set shellmode [Ng_GetCommandLineParameter shellmode]

if { $shellmode == "defined" } {
  set batchmode "defined"
}


if { $batchmode != "defined" } {
    catch {
	wm withdraw .
     
	wm title . $progname
	wm geometry . =850x600
	wm minsize . 400 300
    }
}


source ${ngdir}/variables.tcl
source ${ngdir}/parameters.tcl

if { $batchmode != "defined" } {
    catch {
	source ${ngdir}/menustat.tcl
    }
}

catch { 
    source ${ngdir}/dialog.tcl
}

catch {
    source ${ngdir}/drawing.tcl
}


# if { [catch { load libgeom2dvis[info sharedlibextension] Ng_Geom2d } result ] } {
#    puts "cannot load 2d meshing module" 
#    puts "error: $result"
# }

catch { source ${ngdir}/csgeom.tcl }
catch { source ${ngdir}/stlgeom.tcl }

set hasocc no
catch { source ${ngdir}/occgeom.tcl }

source ${ngdir}/acisgeom.tcl


catch { source ${ngdir}/nghelp.tcl }
catch { source ${ngdir}/ngvisual.tcl }
catch { source ${ngdir}/sockets.tcl }
catch { source ${ngdir}/acis.tcl }



set zugstange 0
catch { source ${ngdir}/trafo/menu.tcl } 



setgranularity ${meshoptions.fineness}

Ng_SetMeshingParameters
Ng_SetVisParameters
Ng_SetDebugParameters
Ng_STLDoctor
Ng_GeometryOptions set

if { $hasocc == "yes" } {
    Ng_SetOCCVisParameters
}


if { $batchmode != "defined" } {
    catch { 
	wm protocol . WM_DELETE_WINDOW { .ngmenu.file invoke "Quit" }
	wm deiconify .
    }
}

set trafoapp 0
catch { source ${ngdir}/trafoapp/trafoapp.tcl }

set geofilename [Ng_GetCommandLineParameter geofile]

if { $geofilename != "undefined" && 
     [info exists trafo] == 0 && $zugstange == 0} {

    if { [ catch { Ng_LoadGeometry $geofilename } errstring] == 0 } {
	if { $batchmode != "defined" } {
	    AddRecentFile $geofilename
	}
	Ng_ParseGeometry
	if { $batchmode != "defined" } {
	    set selectvisual geometry
	    Ng_SetVisParameters
	    redraw
	    wm title . [concat "$progname - " $geofilename]
	}
	set dirname [file dirname $geofilename]
	set basefilename [file tail [file rootname $geofilename]]
    } {
	puts "Problem with input file:"
	puts "$errstring"
    }
}


set cnt 0
foreach { gran } { verycoarse coarse moderate fine veryfine } {
    set cnt [expr $cnt + 1]
    if { [Ng_GetCommandLineParameter $gran] == "defined" } {
	set meshoptions.fineness $cnt
	setgranularity ${meshoptions.fineness}
    }
}


set meshfilename [Ng_GetCommandLineParameter meshfile]
if { $meshfilename == "undefined" } {
    set meshfilename out.mesh
}

set meshfiletype [Ng_GetCommandLineParameter meshfiletype]
if { $meshfiletype == "undefined" } {
    set meshfiletype netgen
}

set inputmeshfilename [Ng_GetCommandLineParameter inputmeshfile]

set mergemeshfilename [Ng_GetCommandLineParameter mergefile]

set meshsizefilename [Ng_GetCommandLineParameter meshsizefile]

if { $meshsizefilename != "undefined" } {
    set options.meshsizefilename $meshsizefilename
}

set refinementfilename [Ng_GetCommandLineParameter refinementfile]


if { $batchmode == "defined" && $solvemode != "defined"} {
    set options.parthread 0
    if { $shellmode == "undefined" } {
# old batchmode: only processes commandline arguments    
      set selectvisual mesh
      Ng_SetVisParameters

      set meshsize [Ng_GetCommandLineParameter meshsize]
      if {$meshsize != "undefined"} { set options.meshsize $meshsize }
        
      if { $inputmeshfilename == "undefined" } {
	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}
      } else {
	Ng_LoadMesh $inputmeshfilename
	if { $mergemeshfilename != "undefined" } {
	    Ng_MergeMesh $mergemeshfilename
        }
      }
	
      if { $refinementfilename != "undefined" } {
	  Ng_Bisect $refinementfilename
      }

      if { $meshfiletype == "netgen" } {
	Ng_SaveMesh $meshfilename
      } else {
	if { [catch { Ng_ExportMesh $meshfilename $meshfiletype } ] == 1 } {
	    puts "Unknown file format $meshfiletype"
        }
      }
      Ng_Exit;

      exit
  } else {
      set code [catch { source ${ngdir}/ngshell.tcl } errcode]
      if {$code} {
	  puts "error: $errcode"
      }  
      set code [ catch {Ng_RunShell} errcode]
      if {$code} {
	  puts "error: $errcode"
      }  
      
      Ng_Exit;
      exit
  }
    
}

set stereo [Ng_GetCommandLineParameter stereo]
if { $stereo == "defined" } {
    set viewoptions.stereo 1 
    puts "use stereo mode" 
    Ng_SetVisParameters; 
    redraw 
}


catch { source ${ngdir}/ngsolve.tcl } 


set scriptfilename [Ng_GetCommandLineParameter script]
if { $scriptfilename != "undefined" } {
    if { [catch { source $scriptfilename } errstring] == 1 } {
	puts "Error in input: $errstring"
    }
}


if { [Ng_GetCommandLineParameter help]=="defined" } {
    if { $zugstange == 1 } {
	print_zug_commandline_help
	exit;
    } {
	if { $trafoapp == 1 } {
	    print_trafo_commandline_help;
	} {
	    print_commandline_help; 
	    Ng_Exit;
	    exit
	}
    }
}


if { [file exists startup.tcl] } {
    source startup.tcl }




catch { source ${ngdir}/demoapp.tcl } 
catch { source ${ngdir}/dropsexp.tcl } 

