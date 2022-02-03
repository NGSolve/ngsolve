puts "loading ngsolve library"

# netgen_library_dir is set from python in pip packages
if { [ info exists netgen_library_dir ] } {
    load $netgen_library_dir/libngsolve[info sharedlibextension] ngsolve
} elseif { [catch { load libngsolve[info sharedlibextension] ngsolve } result ] } {
  set current_script_dir [file dirname [dict get [info frame 0] file]]
  if { [catch { load $current_script_dir/@BIN_TO_LIB_RELPATH@/libngsolve[info sharedlibextension] ngsolve } result2 ] } {
    puts "cannot load ngsolve" 
    puts "error1: $result"
    puts "error2: $result2"
  }
} 

# check some NGS command
if { [catch { NGS_GetData } ] == 0 } { 

    set sbatchmode [Ng_GetCommandLineParameter batchmode]
    if { $sbatchmode == "defined" } {
	puts "batchmode ..."
	
	set pdefilename [Ng_GetCommandLineParameter pdefile]
	if { $pdefilename != "undefined" } {
	    NGS_LoadPDE  $pdefilename;  
	    set solve [Ng_GetCommandLineParameter solve]
	    
	    if { $solve == "defined" } { set solve 1 }
	    if { $solve != "undefined" } {
		set options.parthread 0
		Ng_SetMeshingParameters
		for { set l 1 } { $l <= $solve } { incr l } { NGS_SolvePDE $l }
		NGS_PrintTiming
		puts "Thank you for using $progname/NGSolve"; 
                
                if { [catch { unload libngsolve[info sharedlibextension] ngsolve } result ] } {
                    puts "cannot unload ngsolve" 
                    puts "error: $result"
                } 

		Ng_Exit; 
		destroy .
		exit
	    }
	}
    }
    
    set pdefilename [Ng_GetCommandLineParameter pdefile]
    if { $pdefilename != "undefined" } {
	NGS_LoadPDE  $pdefilename;  

	set solve [Ng_GetCommandLineParameter solve]
	
	if { $solve == "defined" } { set solve 1 }
	if { $solve != "undefined" } {
	    set options.parthread 0
	    Ng_SetMeshingParameters
	    for { set l 1 } { $l <= $solve } { incr l } { NGS_SolvePDE $l }
	    set selectvisual solution
	    Ng_SetVisParameters	
	    redraw
	}
    }
    
    set pyfilename [Ng_GetCommandLineParameter py]
    if { $pyfilename != "undefined" } {
	NGS_LoadPy  $pyfilename;  
    }

    set progname "NGSolve"
    wm title . $progname
    
    .ngmenu add cascade -label "Solve" -menu .ngmenu.solve -underline 1
    
    
    menu .ngmenu.solve
    .ngmenu.solve add command -label "Print Equations" \
	-command { NGS_PrintRegistered }



    menu .ngmenusolvehelp
    .ngmenu.solve add cascade -label "Help" -menu .ngmenusolvehelp 
    
    .ngmenusolvehelp add command -label "Constants..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  constant ] -type ok }
    .ngmenusolvehelp add command -label "Coefficient..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  coefficient ] -type ok }
    .ngmenusolvehelp add command -label "Bilinear-form..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  bilinearform ] -type ok }
    .ngmenusolvehelp add command -label "Linear-form..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  linearform ] -type ok }

    # defined in C++ code:
    menu .ngmenusolvehelpnp
    .ngmenusolvehelp add cascade -label "Numprocs..." -menu .ngmenusolvehelpnp 


    proc SetNumProcHelpMenu {} {
	set allnp [ NGS_Help numprocs ]
	.ngmenusolvehelpnp delete 1 end

	foreach np $allnp {
	    .ngmenusolvehelpnp add command -label "numproc $np" \
		-command "tk_messageBox -title \"Help\" -message \[ NGS_Help numproc $np \] -type ok"
        }
    }

    SetNumProcHelpMenu


    .ngmenu.solve add command -label "Load PDE..." -accelerator "<l><p>"\
	-command { 
	    set types { {"Partial Differential Equation"   {.pde}	} }
	    set file [tk_getOpenFile -filetypes $types -initialdir $dirname]
	    if {$file != ""} {
		set dirname [file dirname $file]
		AddRecentNGSFile $file;
                NGS_LoadPDE  $file;  
		SetNumProcHelpMenu
                set selectvisual mesh;
                Ng_SetVisParameters	
	    }
	}
    
    .ngmenu.solve add cascade -label "Recent Files" -menu .ngmenu.solve.recent 
    menu .ngmenu.solve.recent -tearoff 0

    .ngmenu.solve add command -label "Load Python..." -accelerator "<l><y>"\
	-command { 
	    set types { {"Python script"   {.py}	} }
	    set file [tk_getOpenFile -filetypes $types -initialdir $dirname]
	    if {$file != ""} {
		set dirname [file dirname $file]
		AddRecentPYNGSFile $file;
                NGS_LoadPy  $file;  
	    }
	}


    .ngmenu.solve add command -label "Components..." \
	-command { componentsdialog }
    

    .ngmenu.solve add command -label "Print Report" \
	-command { NGS_PrintPDE }

    .ngmenu.solve add command -label "Memory Usage" \
	-command { NGS_PrintMemoryUsage }

    .ngmenu.solve add command -label "Print Timing" \
	-command { NGS_PrintTiming }


    .ngmenu.solve add command -label "Solve Recent PDE" -accelerator "<s><r>"\
	-command {
#            set filename [.ngmenu.solve.recent entrycget 0 -label]
#            puts "solve filename: $filename"
#	    NGS_LoadPDE  [.ngmenu.solve.recent entrycget 0 -label]
            .ngmenu.solve.recent invoke 0
            if { [ string match "*.pde" [.ngmenu.solve.recent entrycget 0 -label]] } {
                NGS_SolvePDE
            }
	    set selectvisual solution
	    Ng_SetVisParameters	
	    redraw
	}
    
    ttk::button .bubar.pde -text "Recent" \
        -command { .ngmenu.solve invoke "Solve Recent PDE"; }
    pack .bubar.pde -side right
    
    ttk::button .bubar.solve -text "Solve PDE" \
        -command { .ngmenu.solve invoke "Solve PDE"; }
    pack .bubar.solve -side right
    
    ttk::button .bubar.visualize -text "Visual" \
	-command { visual_dialog }
    pack .bubar.visualize -side right
    
    .ngmenu.solve add command -label "Solve PDE" -accelerator "<s><p>"\
	-command {
	    
	    # 	    if {[winfo exists .visoptions_dlg] == 1} {
	    # 		destroy .visoptions_dlg;
	    # 		set test 1;
	    # 	    }
	    
	    
	    NGS_SolvePDE
	    set selectvisual solution
	    Ng_SetVisParameters	

	    Ng_Vis_Set parameters; 

	    #	    if { $test == 1 } {
	    #		reset_visual_dialog 
	    #	    }
	    
	    redraw
	}

    .ngmenu.solve add cascade -label "Solve PDE x" -menu .ngmenu.solve.solvex
    menu .ngmenu.solve.solvex

    proc SolveX { num } {
	for { set i 1 } { $i <= $num } { incr i } {
	    uplevel 1 "NGS_SolvePDE $i"
	}
    }

    .ngmenu.solve.solvex  add command -label "1 Level" -command { SolveX 1 }
    .ngmenu.solve.solvex  add command -label "2 Level" -command { SolveX 2 }
    .ngmenu.solve.solvex  add command -label "3 Level" -command { SolveX 3 }
    .ngmenu.solve.solvex  add command -label "4 Level" -command { SolveX 4 }
    .ngmenu.solve.solvex  add command -label "5 Level" -command { SolveX 5 }
    .ngmenu.solve.solvex  add command -label "6 Level" -command { SolveX 6 }
    .ngmenu.solve.solvex  add command -label "7 Level" -command { SolveX 7 }
    .ngmenu.solve.solvex  add command -label "8 Level" -command { SolveX 8 }
    .ngmenu.solve.solvex  add command -label "9 Level" -command { SolveX 9 }
    .ngmenu.solve.solvex  add command -label "10 Level" -command { SolveX 10 }
    .ngmenu.solve.solvex  add command -label "11 Level" -command { SolveX 11 }
    .ngmenu.solve.solvex  add command -label "12 Level" -command { SolveX 12 }
    .ngmenu.solve.solvex  add command -label "13 Level" -command { SolveX 13 }
    .ngmenu.solve.solvex  add command -label "14 Level" -command { SolveX 14 }
    .ngmenu.solve.solvex  add command -label "15 Level" -command { SolveX 15 }
    .ngmenu.solve.solvex  add command -label "16 Level" -command { SolveX 16 }
    .ngmenu.solve.solvex  add command -label "17 Level" -command { SolveX 17 }
    .ngmenu.solve.solvex  add command -label "18 Level" -command { SolveX 18 }
    .ngmenu.solve.solvex  add command -label "19 Level" -command { SolveX 19 }
    .ngmenu.solve.solvex  add command -label "20 Level" -command { SolveX 20 }


    .ngmenu.solve add command -label "Enter Command" \
	-command { 
	    NGS_EnterCommand;
	}


    .ngmenu.solve add command -label "Visualization..." \
	-command { 
	    visual_dialog;
	}
    
    
    .ngmenu.solve add command -label "Save Solution..." \
	-command { 
	    set types { {"Solution File"  {.sol} } }
	    set file [tk_getSaveFile -filetypes $types -defaultextension ".sol"  ]
	    if {$file != ""} {
		NGS_SaveSolution $file 
	    }
	}
    
    .ngmenu.solve add command -label "Load Solution..." \
	-command { 
	    set types { {"Solution File"  {.sol} } }
	    set file [tk_getOpenFile -filetypes $types -defaultextension ".sol"  ]
	    if {$file != ""} {
		NGS_LoadSolution $file 
		set selectvisual solution
		Ng_SetVisParameters
		redraw
	    }
	}

    .ngmenu.solve add command -label "Dump PDE" \
        -command {
            NGS_DumpPDE test.arch
        }

    .ngmenu.solve add command -label "Restore PDE" \
        -command {
            NGS_RestorePDE test.arch
        }

    .ngmenu.solve add command -label "Socket-load" \
	-command { 
            NGS_SocketLoad 52002  128.131.37.12
            # NGS_SocketLoad 52002  localhost
            # numericus.asc.tuwien.ac.at
            set selectvisual solution
	    Ng_SetVisParameters	
	    redraw
        }

    .ngmenu.solve add command -label "Python shell"  -accelerator "<p><y>"\
	-command { 
            NGS_PythonShell
	}



    #     .ngmenu.solve add command -label "Save Solution (ASCII)..." \
	# 	-command { 
    # 	    set types { {"Solution File"  {.asol} } }
    # 	    set file [tk_getSaveFile -filetypes $types -defaultextension ".asol"  ]
    # 	    if {$file != ""} {
    # 		NGS_SaveSolution $file 1
    # 	    }
    # 	}
    
    #     .ngmenu.solve add command -label "Load Solution (ASCII)..." \
	# 	-command { 
    # 	    set types { {"Solution File"  {.asol} } }
    # 	    set file [tk_getOpenFile -filetypes $types -defaultextension ".asol"  ]
    # 	    if {$file != ""} {
    # 		NGS_LoadSolution $file 1
    # 		set selectvisual solution
    # 		Ng_SetVisParameters
    # 		redraw
    # 	    }
    # 	}



    #    .ngmenu.solve add separator



    .ngmenu.help delete "About..."
    .ngmenu.help add command -label "About..." \
	-command {
	    tk_messageBox -message \
		"This is NETGEN/NGSolve \n mainly written by \n Joachim Schoeberl \n\
                 at Vienna University of Technology, Austria \n\
                 and RWTH Aachen University, Germany \n\
                 and Johannes Kepler University, Linz, Austria \n\
                 supported by the Austrian Science Foundation FWF \n\
                 thanks to \n\
                 F. Bachinger, A. Becirovic, H. Egger, R. Gaisbauer, J. Gerstmayr, M. Hochsteger, G. Kitzler, U. Langer, C. Lehrenfeld, P. Rajan, A. Sinwel, M. Wabro, S. Zaglmayr"
	}

    proc AddRecentNGSFile { filename } {
	global progname
	catch { [.ngmenu.solve.recent delete $filename] }
	.ngmenu.solve.recent insert 0 command -label $filename \
	    -command "AddRecentNGSFile {$filename}; 
		NGS_LoadPDE  {$filename};
                set selectvisual mesh;
                Ng_SetVisParameters	
                wm title . [concat \" $progname - $filename \"];"
	
	if { [.ngmenu.solve.recent index last] >= 6 } {
	    .ngmenu.solve.recent delete last }

	savengsinifile;
    }
    proc AddRecentPYNGSFile { filename } {
	global progname
	catch { [.ngmenu.solve.recent delete $filename] }
	.ngmenu.solve.recent insert 0 command -label $filename \
	    -command "AddRecentPYNGSFile {$filename}; 
                wm title . [concat \" $progname - $filename \"];
		NGS_LoadPy  {$filename};"
	
	if { [.ngmenu.solve.recent index last] >= 6 } {
	    .ngmenu.solve.recent delete last }

	savengsinifile;
    }
    

    # the ini file is saved on demand :
    global ngsinifilename

    set ngsinifilename [file join $nguserdir ngs.ini]
    
    proc savengsinifile { } {
	global ngsinifilename
	if {[catch { set datei [open $ngsinifilename w]  } result ]} {
	    puts "cannot write to $ngsinifilename file"
	} {
	    for { set i [.ngmenu.solve.recent index last] } { $i >= 0 } { incr i -1 } {
		puts $datei "recentfile \"[.ngmenu.solve.recent entrycget $i -label]\""
	    }
	    close $datei
	}
    }
    
    proc loadngsinifile { } {
	global ngsinifilename
	if { [file exists $ngsinifilename] == 1 } {
	    set datei [open $ngsinifilename r]
	    while { [gets $datei line] >= 0 } {
		if {[lindex $line 0] == "recentfile"} {
                    if { [string match "*.pde" [lindex $line 1]] } {
                        AddRecentNGSFile [lindex $line 1];
                    }
                    if { [string match "*.py" [lindex $line 1]] } {
                        AddRecentPYNGSFile [lindex $line 1];
                    }
		}
	    }
	    close $datei
	}
    }

    loadngsinifile;


    proc componentsdialog { } {
	
	set w .components_dlg
	
	if {[winfo exists .components_dlg] == 1} {
	    wm withdraw $w
	    wm deiconify $w
	    focus $w 
	} {

	    toplevel $w
        #tk::treeview $w.tree
        #pack $w.tree -fill both -expand y
        #$w.tree insert {} end -id widgets -text "Widget Tour"
        ttk::treeview $w.tree
        pack $w.tree -fill both -expand y
        $w.tree insert {} end -id constants -text "Constants"
        $w.tree insert {} end -id variables -text "Variables"
        $w.tree insert {} end -id coefficients -text "Coefficients"
        $w.tree insert {} end -id spaces -text "Spaces"
        $w.tree insert {} end -id biforms -text "Bilinear-forms"
        $w.tree insert {} end -id liforms -text "Linear-forms"
        $w.tree insert {} end -id gfs -text "Grid-functions"
        $w.tree insert {} end -id preconds -text "Preconditioners"
        $w.tree insert {} end -id nps -text "NumProcs"
        set constants [NGS_GetData constants]
	    foreach co $constants {
        $w.tree insert {constants} end -text $co
	    }
        set variables [NGS_GetData variableswithval]
	    foreach va $variables {
        $w.tree insert {variables} end -text $va
	    }
	    set coefs [NGS_GetData coefficients]
	    foreach coef $coefs {
		$w.tree insert {coefficients} end -text $coef
	    }        
	    set spaces [NGS_GetData spaces]
	    foreach space $spaces {
		$w.tree insert {spaces} end -text $space
	    }
	    set biforms [NGS_GetData bilinearforms]
	    foreach biform $biforms {
		$w.tree insert {biforms} end -text $biform
	    }
	    set liforms [NGS_GetData linearforms]
	    foreach liform $liforms {
		$w.tree insert {liforms} end -text $liform
	    }
	    set gridfuns [NGS_GetData gridfunctions]
	    foreach gridfun $gridfuns {
		$w.tree insert {gfs} end -text $gridfun
	    }
	    set preconds [NGS_GetData preconditioners]
	    foreach precond $preconds {
		$w.tree insert {preconds} end -text $precond
	    }
	    set numprocs [NGS_GetData numprocs]
	    foreach numproc $numprocs {
		$w.tree insert {nps} end -text $numproc        
        }
        #$w.tree configure -color black
	    # tixTree $w.mtre -options { separator "\\" }
	    # pack $w.mtre -fill both -expand y
	    # set hlist [$w.mtre subwidget hlist]

            # $hlist configure -selectforeground black
            # $hlist configure -selectbackground grey

	    # $hlist add constants -itemtype text -text "Constants"
	    # set constants [NGS_GetData constants]
	    # foreach co $constants {
		# $hlist add constants\\$co -itemtype text -text $co
	    # }

	    # $hlist add variables -itemtype text -text "Variables"
	    # set variables [NGS_GetData variableswithval]
	    # foreach va $variables {
		# $hlist add variables\\$va -itemtype text -text $va
	    # }

	    # $hlist add coeffs -itemtype text -text "Coefficients"
	    # set coefs [NGS_GetData coefficients]
	    # foreach coef $coefs {
		# $hlist add coeffs\\$coef -itemtype text -text $coef
	    # }


	    # $hlist add spaces -itemtype text -text "Spaces"
	    # set spaces [NGS_GetData spaces]
	    # foreach space $spaces {
		# $hlist add spaces\\$space -itemtype text -text $space
	    # }

	    # $hlist add biforms -itemtype text -text "Bilinear-forms"
	    # set biforms [NGS_GetData bilinearforms]
	    # foreach biform $biforms {
		# $hlist add biforms\\$biform -itemtype text -text $biform
	    # }

	    # $hlist add liforms -itemtype text -text "Linear-forms"
	    # set liforms [NGS_GetData linearforms]
	    # foreach liform $liforms {
		# $hlist add liforms\\$liform -itemtype text -text $liform
	    # }

	    # $hlist add gridfuns -itemtype text -text "Grid-functions"
	    # set gridfuns [NGS_GetData gridfunctions]
	    # foreach gridfun $gridfuns {
		# $hlist add gridfuns\\$gridfun -itemtype text -text $gridfun
	    # }

	    # $hlist add preconds -itemtype text -text "Preconditioners"
	    # set preconds [NGS_GetData preconditioners]
	    # foreach precond $preconds {
		# $hlist add preconds\\$precond -itemtype text -text $precond
	    # }

	    # $hlist add numprocs -itemtype text -text "NumProcs"
	    # set numprocs [NGS_GetData numprocs]
	    # foreach numproc $numprocs {
		# $hlist add numprocs\\$numproc -itemtype text -text $numproc
	    # }

	    #	    $hlist add varis -itemtype text -text "Variables"
	    #	    set varis [NGS_GetData variables]
	    #	    foreach vari $varis {
	    #		$hlist add varis\\$vari -itemtype text -text $vari
	    #	    }

	    # $hlist add varis2 -itemtype text -text "Variables2"
	    # set varis2 [NGS_GetData variablesval]
	    # foreach vari $varis2 {
	    # scan $vari "val%fname%s" vval vname
	    # $hlist add varis2\\$vari  -itemtype text -text [format "%s %1.3e" $vname $vval]
	    # }

	    



	    # 	    for { set i 1 } { $i <= $nspaces } { incr i } {
	    # 		set name [NGS_GetData spacename $i]
	    # 		$hlist add spaces\\$name -itemtype text -text $name
	    # 	    }

	    # 	    $hlist add gridfunctions -itemtype text -text "Grid-functions"
	    # 	    set ngf [NGS_GetData numgridfunctions]
	    # 	    for { set i 1 } { $i <= $ngf } { incr i } {
	    # 		set name [NGS_GetData gridfunctionname $i]
	    # 		$hlist add gridfunctions\\$name -itemtype text -text $name
	    # 	    }

	    # 	    $hlist add bf -itemtype text -text "Bilinear-forms"
	    # 	    set nbf [NGS_GetData numbilinearforms]
	    # 	    for { set i 1 } { $i <= $nbf } { incr i } {
	    # 		set name [NGS_GetData bilinearformname $i]
	    # 		$hlist add bf\\$name -itemtype text -text $name

	    # 		set nbfi [NGS_GetData numbilinearformcomps $name]
	    # 		for { set j 1 } { $j <= $nbfi } { incr j } {
	    # 		    set compname [NGS_GetData bilinearformcompname $name $j]
	    # 		    $hlist add bf\\$name\\$j -itemtype text -text $compname
	    # 		}
	    # 	    }

	    # 	    $hlist add lf -itemtype text -text "Linear-forms"
	    # 	    set nlf [NGS_GetData numlinearforms]
	    # 	    for { set i 1 } { $i <= $nlf } { incr i } {
	    # 		set name [NGS_GetData linearformname $i]
	    # 		$hlist add lf\\$name -itemtype text -text $name

	    # 		set nlfi [NGS_GetData numlinearformcomps $name]
	    # 		for { set j 1 } { $j <= $nlfi } { incr j } {
	    # 		    set compname [NGS_GetData linearformcompname $name $j]
	    # 		    $hlist add lf\\$name\\$j -itemtype text -text $compname
	    # 		}
	    # 	    }

        # Old version
	    # $w.mtre autosetmode
	    

	    # bind $hlist <Double-1> {
		# set solname [[.components_dlg.mtre subwidget hlist] info selection]
		# puts $solname
		# set seppos [string first \\ $solname]
		# if { $seppos != -1 } {
		    # set field [string range $solname 1 [expr $seppos-1]]
		    # set name [string range $solname [expr $seppos+1] [expr [string length $solname]-2]]
		    # puts "field = $field, name = $name"
		    # NGS_PrintPDE $field $name
		# }
	    # }

            # old version end
            #new version
            bind $w.tree <Double-1> {
                set name [.components_dlg.tree item [.components_dlg.tree selection] -text ]
                # set field [.components_dlg.tree item [.components_dlg.tree parent [.components_dlg.tree selection] ] -text]
                set field [.components_dlg.tree parent [.components_dlg.tree selection] ]
                
                NGS_PrintPDE $field $name
                # puts "field = $field, name = $name"
                
            }
            #new version end
        
	    button $w.cl -text "Close" -command {
		destroy .components_dlg
	    }

	    pack  $w.cl
	    
	    
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "Components"
	    focus .components_dlg
	}
    }




    if { [Ng_GetCommandLineParameter recent]=="defined" } {
	if { [catch { .ngmenu.solve invoke "Solve Recent PDE";  } errstring] == 1 } {
	    puts "TCL-ERROR handler:\n $errstring";
	    exit;
	}
    }






    bind . <l><p> { .ngmenu.solve invoke "Load PDE..." }  ; 
    bind . <l><y> { .ngmenu.solve invoke "Load Python..." }  ; 
    bind . <s><r> { .ngmenu.solve invoke "Solve Recent PDE" }  ; 
    bind . <s><p> { .ngmenu.solve invoke "Solve PDE" }  ; 
    bind . <e><c> { .ngmenu.solve invoke "Enter Command" }  ; 
    bind . <p><y> { .ngmenu.solve invoke "Python shell" }  ; 


}

