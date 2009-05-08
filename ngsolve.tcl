# check some NGS command
if { [catch { NGS_GetData } ] == 0 } { 
    
    set progname "NGSolve"
    wm title . $progname
    
    .ngmenu add cascade -label "Solve" -menu .ngmenu.solve -underline 1
    
    # menubutton .mbar.solve -text Solve -underline 0 -menu .mbar.solve.menu
    # pack .mbar.solve -side left
    
    menu .ngmenu.solve
    .ngmenu.solve add command -label "Print Equations" \
	-command { NGS_PrintRegistered }





    menu .ngmenusolvehelp
    .ngmenu.solve add cascade -label "Help" -menu .ngmenusolvehelp
    
    .ngmenusolvehelp add command -label "Coefficient..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  coefficient ] -type ok }
    .ngmenusolvehelp add command -label "Bilinear-form..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  bilinearform ] -type ok }
    .ngmenusolvehelp add command -label "Linear-form..." \
	-command { tk_messageBox -title "Help" -message  [ NGS_Help  linearform ] -type ok }

    # defined in C++ code:
    .ngmenusolvehelp add cascade -label "Numprocs..." -menu .ngmenusolvehelpnp 
    
    .ngmenusolvehelp add command -label "Latest News..."  \
	-command { tk_messageBox -title "Latest News" -message \
                       { 
                           06042004 online documentation (JS) 
                       } -type ok }  ;
    


    .ngmenu.solve add command -label "Load PDE..." -accelerator "<l><p>"\
	-command { 
	    set types { {"Partial Differential Equation"   {.pde}	} }
	    set file [tk_getOpenFile -filetypes $types]
	    if {$file != ""} {
		AddRecentNGSFile $file;
                NGS_LoadPDE  $file;  
                set selectvisual mesh;
                Ng_SetVisParameters	
	    }
	}
    
    .ngmenu.solve add cascade -label "Recent Files" -menu .ngmenu.solve.recent 
    menu .ngmenu.solve.recent

    .ngmenu.solve add command -label "Components..." \
	-command { componentsdialog }
    

    .ngmenu.solve add command -label "Print Report" \
	-command { NGS_PrintPDE }

    .ngmenu.solve add command -label "Memory Usage" \
	-command { NGS_PrintMemoryUsage }

    .ngmenu.solve add command -label "Print Timing" \
	-command { NGS_PrintTiming }


#    .ngmenu.solve add command -label "Play anim" \
# 	-command { 
# 	    set selectvisual mesh;
# 	    for { set i 1 } { $i <= 1000 } { incr i } {
# 		#	    NGS_PlayAnim 5 175 kurbc;
# 		#	    NGS_PlayAnim kurb 2 500;
# 		NGS_PlayAnim 1 50 kurbc 
# 		update 
# 		.ndraw render;
# 	    }

	    #	set types { {"MBS Solution Data"   {.sol}	} }
	    #	set file [tk_getOpenFile -filetypes $types]
	    #	if {$file != ""} {
	    #	    NGS_PlayAnim  $file;
	    #	}
#	}




    .ngmenu.solve add command -label "Solve Recent PDE" -accelerator "<s><r>"\
	-command { 
	    NGS_LoadPDE  [.ngmenu.solve.recent entrycget 1 -label]
	    NGS_SolvePDE
	    set selectvisual solution
	    Ng_SetVisParameters	


	    redraw
	}

#    .ngmenu.solve add command -label "Test Reissner-Mindlin" \
\#	-command { 
#	    NGS_TestRM
#	}



    button .bubar.pde -text "Recent" \
	-command { .ngmenu.solve invoke "Solve Recent PDE"; }
    pack .bubar.pde -side right

    button .bubar.solve -text "Solve" \
	-command { .ngmenu.solve invoke "Solve PDE"; }
    pack .bubar.solve -side right

    button .bubar.visualize -text "Visual" \
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



#    for { set lev 1 } { $lev <= 20 } { incr lev } {
#	set a { puts $lev; }
#	.ngmenu.solve.solvex  add command -label "$lev Level" \
\#	    -command { $a }
#    }
#		for { set i 1 } { $i <= $lev } { incr i } {
#		    puts "i = $i, lev = [expr $lev]"
#		    NGS_SolvePDE
#    }


    #    .ngmenu.solve add command -label "Solve MBS" 
    #  -command  { NGS_MBSSolver "examples/kurbel.mbs" }
    #	set types { {"MBS File"   {.mbs}	} }
    #	set file [tk_getOpenFile -filetypes $types]
    #	if {$file != " "} {
    #	    NGS_MBSSolver $file;  }


    .ngmenu.solve add command -label "Visualization..." \
	-command { 
	    visual_dialog;
	}
    

#     .ngmenu.solve add command -label "Demo session" \
# 	-command {
# 	    tk_messageBox -message "Load PDE file";
# 	    NGS_LoadPDE ngsolve/pde_tutorial/d7_coil.pde;

# 	    tk_messageBox -message "Solve PDE";
# 	    NGS_SolvePDE

# 	    tk_messageBox -message "Switch to solution";
# 	    set selectvisual solution
# 	    Ng_SetVisParameters	
# 	    redraw

# 	    tk_messageBox -message "Define clipping plane";
# #	    clippingdialog
# 	    set viewoptions.clipping.enable 1
# 	    Ng_SetVisParameters; 
# 	    redraw;

# 	    tk_messageBox -message "Open visualization dialog"
# 	    visual_dialog;
# 	    tk_messageBox -message "Define clipping plane";
# 	}


    
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
            "This is NETGEN/NGSolve \n mainly written by \n Joachim Schöberl \n\
                 at RWTH Aachen University, Germany \n\
                 and Johannes Kepler University, Linz, Austria \n\
                 supported by the Austrian Science Foundation FWF \n\
                 thanks to \n\
                 F. Bachinger, A. Becirovic, H. Egger, R. Gaisbauer, J. Gerstmayr, U. Langer, A. Sinwel, M. Wabro, S. Zaglmayr"
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
    

    # the ini file is saved  on demand :
    proc savengsinifile { } {
	uplevel 1  {
            if {[catch { set datei [open ngs.ini w]  } result ]} {
                puts "cannot write to ng.ini file"
            } {
                for { set i [.ngmenu.solve.recent index last] } { $i >= 1 } { incr i -1 } {
                    puts $datei "recentfile \"[.ngmenu.solve.recent entrycget $i -label]\""
                }
                close $datei
            }
	}
    }
    
    proc loadngsinifile { } {
	if { [file exists ngs.ini] == 1 } {
	    set datei [open ngs.ini r]
	    while { [gets $datei line] >= 0 } {
		if {[lindex $line 0] == "recentfile"} {
		    AddRecentNGSFile [lindex $line 1]
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


	    tixTree $w.mtre -options { separator "\\" }
	    pack $w.mtre -fill both -expand y
	    set hlist [$w.mtre subwidget hlist]




	    $hlist add coeffs -itemtype text -text "Coefficients"
	    set coefs [NGS_GetData coefficients]
	    foreach coef $coefs {
		$hlist add coeffs\\$coef -itemtype text -text $coef
	    }


 	    $hlist add spaces -itemtype text -text "Spaces"
	    set spaces [NGS_GetData spaces]
	    foreach space $spaces {
		$hlist add spaces\\$space -itemtype text -text $space
	    }

 	    $hlist add biforms -itemtype text -text "Bilinear-forms"
	    set biforms [NGS_GetData bilinearforms]
	    foreach biform $biforms {
		$hlist add biforms\\$biform -itemtype text -text $biform
	    }

 	    $hlist add liforms -itemtype text -text "Linear-forms"
	    set liforms [NGS_GetData linearforms]
	    foreach liform $liforms {
		$hlist add liforms\\$liform -itemtype text -text $liform
	    }

 	    $hlist add gridfuns -itemtype text -text "Grid-functions"
	    set gridfuns [NGS_GetData gridfunctions]
	    foreach gridfun $gridfuns {
		$hlist add gridfuns\\$gridfun -itemtype text -text $gridfun
	    }

 	    $hlist add preconds -itemtype text -text "Preconditioners"
	    set preconds [NGS_GetData preconditioners]
	    foreach precond $preconds {
		$hlist add preconds\\$precond -itemtype text -text $precond
	    }

 	    $hlist add numprocs -itemtype text -text "NumProcs"
	    set numprocs [NGS_GetData numprocs]
	    foreach numproc $numprocs {
		$hlist add numprocs\\$numproc -itemtype text -text $numproc
	    }

#	    $hlist add varis -itemtype text -text "Variables"
#	    set varis [NGS_GetData variables]
#	    foreach vari $varis {
#		$hlist add varis\\$vari -itemtype text -text $vari
#	    }

#	    $hlist add varis2 -itemtype text -text "Variables2"
#	    set varis2 [NGS_GetData variablesval]
#	    foreach vari $varis2 {
#		scan $vari "val%fname%s" vval vname
#		$hlist add varis2\\$vari  -itemtype text -text [format "%s %1.3e" $vname $vval]
#	    }

	    



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

	    $w.mtre autosetmode
	    

	    bind $hlist <Double-1> {
		set solname [[.components_dlg.mtre subwidget hlist] info selection]
		puts $solname
		set seppos [string first \\ $solname]
		if { $seppos != -1 } {
		    set field [string range $solname 1 [expr $seppos-1]]
		    set name [string range $solname [expr $seppos+1] [expr [string length $solname]-2]]
		    puts "field = $field, name = $name"
		    NGS_PrintPDE $field $name
		}
	    }

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



set pdefilename [Ng_GetCommandLineParameter pdefile]
if { $pdefilename != "undefined" } {
    NGS_LoadPDE  $pdefilename;  

    set solve [Ng_GetCommandLineParameter solve]
    if { $zugstange == 1 } {
	set options.parthread 0
	NGS_SolvePDE;
    } {
	if { $solve == "defined" } {
	    set options.parthread 0
	    NGS_SolvePDE
	    exit;
	} {
	    if { $solve != "undefined" } {
		set options.parthread 0
		for { set l 1 } { $l <= $solve } { incr l } { NGS_SolvePDE $l }
		exit;
	    }
	}
    }
}



bind . <l><p> { .ngmenu.solve invoke "Load PDE..." }  ; 
bind . <s><r> { .ngmenu.solve invoke "Solve Recent PDE" }  ; 
bind . <s><p> { .ngmenu.solve invoke "Solve PDE" }  ; 





}

