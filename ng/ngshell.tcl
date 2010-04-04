## the shell loop

proc dotest {} {
  source ngtest.tcl
}

proc Ng_RunShell {} {
  puts "Wellcome to NG Shell mode"
  set line 1
  while { 1 } {
    puts -nonewline "$line: "
    flush stdout
    set cmdline [gets stdin] 
    if { [catch $cmdline errcode] } {
#      puts "error in command: '$cmdline'"
      puts "$errcode"
    } 
    incr line 1
  }
}

## global list for help index
# ---> global var in variables.tcl
# set cmdindex {}
# set hlpindex {}
# set secindex {}



# print comd list
proc Ng_PrintCmdIndex { } {
  global cmdindex
  foreach { lst } $cmdindex {
    puts $lst 
  }  
}

# print formatted help index
proc Ng_PrintHlpIndex { } {
  global hlpindex
  global secindex
  foreach {sec} $secindex {
    puts "\n  * $sec:"
    foreach {lst} $hlpindex {
      if {$sec == [lindex $lst 1]} {
        puts "    * [lindex $lst 2]: [lindex $lst 3]"
      }
    }
  } 
}

# register a cmd to the help index
proc Ng_RegisterCmd { cmd section syntax {help ""} } {
  global hlpindex
  global cmdindex
  global secindex

    puts "register command $cmd"

  if { [lsearch $cmdindex cmd] != -1 } {
    puts "command '$cmd' already defined"
  } else {
    lappend cmdindex $cmd

    lappend hlpindex [list $cmd $section $syntax $help]
    if {[lsearch $secindex $section]==-1} {
      lappend secindex $section
    }
#    puts "registered command $cmd"
  }
}


# general purpose commands
Ng_RegisterCmd "exit" "general" "exit" "exit Netgen shell mode" 
#Ng_RegisterCmd "Ng_LoadGeometry" "netgen" "Ng_LoadGeometry <file>" "load geometry file" 
#Ng_RegisterCmd "Ng_ParseGeometry" "netgen" "Ng_ParseGeometry" "parse geometry"
#Ng_RegisterCmd "Ng_GenerateMesh" "netgen" "Ng_GenerateMesh" "generate mesh"
#Ng_RegisterCmd "Ng_SaveMesh" "netgen" "Ng_SaveMesh <file>" "save mesh to file" 


## public domain shell functions

# print hel information
proc nghelp { {sec ""} } {
  global secindex
  global hlpindex
  global cmdindex
  if { $sec == "" } {
      Ng_PrintHlpIndex
      puts "\n  type help 'section'\n"
      return
  } 
  
  
  if { [lsearch $secindex $sec] != -1} {
    foreach {lst} $hlpindex {
      if {[lindex $lst 1] == $sec } {
        puts "  * [lindex $lst 2]: [lindex $lst 3]"
      }
    }
    return
  }
  
  set ind [lsearch $cmdindex $sec]
  if {$ind != -1} {
    set lst [lindex $hlpindex $ind]
    puts "  * [lindex $lst 2]: [lindex $lst 3]"
    return
  }
  
  puts "  unknown section or command $sec"
}

set ngtimer 0

proc nggettimer {} {
   return [clock clicks -milliseconds]
}

proc ngtic {} {
  set ::ngtimer [nggettimer]
}


proc ngtoc { {logfile stdout} } {
  set end [nggettimer]
  set tim [expr ($end - $::ngtimer)/1000.0]
  puts $logfile "$tim s"
}


# load geometry file
proc ngloadgeometry { fname } {
  if { ![file exists $fname] } {
    puts "error: file $fname does not exist"
  } else {
    set err [catch {Ng_LoadGeometry $fname}]
    if {$err != 0} {
      puts "error: loading geometry failed"
    }
  }
} 
Ng_RegisterCmd "ngloadgeometry" "netgen" "ngloadgeometry <file>" "load geometry file"  

# parse geometry
proc ngparsegeometry {} {
  set err [catch {Ng_ParseGeometry}]
  if {$err} {
    puts "error: parsing geometry failed"
  }
}
Ng_RegisterCmd "ngparsegeometry" "netgen" "ngparsegeometry" "parse geometry"  
 
# generate mesh
proc nggeneratemesh {} {
  set err [catch {Ng_GenerateMesh}]
  if {$err} {
    puts "error:  mesh generation failed"
  }
}
Ng_RegisterCmd "nggeneratemesh" "netgen" "nggeneratemesh" "generate mesh"  

# save mesh
proc ngsavemesh { fname } {
  if { [file exists $fname]} {
    puts "warning: existing file $fname overwritten"
  } else {
    set err [catch {Ng_SaveMesh $fname}]
    if {$err != 0} {
      puts "error: saving mesh failed"
    }
  }
}  
Ng_RegisterCmd "ngsavemesh" "netgen" "ngsavemesh <file>" "save mesh to file"  

#set option
proc ngset { opt {val 0} } {
  if {$opt == "meshsize"} {
    set ::options.meshsize $val
    Ng_SetMeshingParameters
  } elseif {$opt == "printmsg"} {
    set ::options.printmsg $val
    Ng_SetMeshingParameters
  } else {
    puts "error: unknown option $opt";
  }
}
Ng_RegisterCmd "ngset" "netgen" "ngset <option> <val>" "set option to val"  

proc nganalyzegeometry {} {
  Ng_GenerateMesh ag ag
  Ng_ReadStatus
}

proc ngmeshedges {} {
  Ng_GenerateMesh me me 
  Ng_ReadStatus
}

proc ngmeshsurface { } {
 Ng_GenerateMesh ms ms
 Ng_ReadStatus
}

proc ngoptimizesurface { {step all} } {
  if {$step == "all"} {
    Ng_GenerateMesh os os cmsmSm
  } elseif {$step == "meshsmoothing"} {
    Ng_GenerateMesh os os m
  } elseif {$step == "topologicedgeswapping" } {
    Ng_GenerateMesh os os s
  } elseif {$step == "metricedgeswapping"} {
    Ng_GenerateMesh os os S
  } elseif {$step == "combinepoints"} {
    Ng_GenerateMesh os os c 
  } else {
    puts "error: unkown option in ngoptimizesurface"
    return
  }
  Ng_ReadStatus
}

proc ngmeshvolume { } {
  Ng_GenerateMesh mv mv 
  Ng_ReadStatus 
}

proc ngoptimizevolume {{step ""} } {
  if {$step == ""} {
    Ng_GenerateMesh ov ov 
  } elseif {$step == "smooth"} {
    Ng_GenerateMesh ov ov m
  } elseif {$step == "smoothjacobian"} {
    Ng_GenerateMesh ov ov j 
  } else {
    puts "error: unknown step $step in ngoptimizevolume"
    return
  }
  Ng_ReadStatus 
}


proc ngsloadpde {fname} {
  if { ![file exists $fname] } {
    puts "warning: pdefile $fname does not exist"
  } else {
    puts "load pde $fname"
    NGS_LoadPDE $fname
  }
}

proc ngssolvepde {} { 
  NGS_SolvePDE 
}

catch {source "${::ngdir}/ngtesting.tcl"} errcode
# puts "errcode = $errcode"
