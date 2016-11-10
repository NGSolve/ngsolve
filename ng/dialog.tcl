proc meshingoptionsdialog { } {

    set w .options_dlg
    
    if {[winfo exists .options_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w
        
    } {
        
	toplevel $w
        #wm resizable $w 0 0 

#	global options.meshsize

        pack [ttk::notebook $w.nb]  -fill both -side top
        $w.nb add [ttk::frame $w.nb.general] -text "General" -underline 0 
        $w.nb add [ttk::frame $w.nb.meshsize] -text "Mesh Size" -underline 0
        $w.nb add [ttk::frame $w.nb.chartopt] -text "STL Charts" -underline 0
        $w.nb add [ttk::frame $w.nb.optimizer] -text "Optimizer" -underline 0
        # $w.nb add [ttk::frame $w.nb.insider] -text "Insider" -underline 0
        $w.nb add [ttk::frame $w.nb.debug] -text "Debug" -underline 0
        
	# tixNoteBook $w.nbold -ipadx 6 -ipady 6
	# $w.nbold add general -label "General" -underline 0 
	# $w.nbold add meshsize -label "Mesh Size"   -underline 0
	# $w.nbold add chartopt -label "STL Charts" -underline 0
	# $w.nbold add optimizer -label "Optimizer"   -underline 0
	# $w.nbold add insider -label "Insider"   -underline 0
	# $w.nbold add debug -label "Debug"   -underline 0
	# pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	


        # ############################################################
	# General meshing options
        # ############################################################
        
        set f $w.nb.general
        ttk::frame $f.background
        pack $f.background -fill both 
        set f $f.background
        ttk::labelframe $f.f2 -relief groove -borderwidth 3 -text "General meshing options"
        pack $f.f2  -pady 15 -fill x 
        set f $f.f2
        
	set finevals { 1 2 3 4 5 6 }
	set finelabs(1) "very coarse" 
	set finelabs(2) "coarse" 
	set finelabs(3) "moderate" 
	set finelabs(4) "fine" 
	set finelabs(5) "very fine" 
	set finelabs(6) "user defined"

	#tixOptionMenu $f.fine -label "Mesh granularity : " \
	    -options {
	#	label.width  19
	#	label.anchor e
	#	menubutton.width 15
	#    } 

	#foreach finev $finevals {
	#    $f.fine add command $finev -label $finelabs($finev)
	#}
	#$f.fine config -variable meshoptions.fineness
	#$f.fine config -command { setgranularity }
	#global meshoptions.fineness
        #setgranularity ${meshoptions.fineness}
        #pack $f.fine

        global meshoptions.fineness
        ttk::label  $f.fine2l -text "Mesh granularity: "
        ttk::menubutton $f.fine2c -menu $f.fine2m -text "coarse" -width 16

        menu $f.fine2m  -tearoff 0
	foreach finev { 1 2 3 4 5 6 } {
	    $f.fine2m add command -label $finelabs($finev) \
                -command "set meshoptions.fineness $finev ; setgranularity $finev; $f.fine2c configure -text \"$finelabs($finev)\""
	}
        $f.fine2m invoke $finelabs(${meshoptions.fineness})                


        grid $f.fine2l $f.fine2c -sticky nw
        

	set mgsteps { ag me ms os mv ov }
	set mgsteplabel(ag) "Analyze Geometry"
	set mgsteplabel(me) "Mesh Edges"
	set mgsteplabel(ms) "Mesh Surface"
	set mgsteplabel(os) "Optimize Surface"
	set mgsteplabel(mv) "Mesh Volume"
	set mgsteplabel(ov) "Optimize Volume"

        global meshoptions.firststep 
        ttk::label  $f.first2l -text "First Step: "
        # ttk::menubutton $f.first2.c -menu $f.first2.m -text "Analyze Geometry" -width 12
        ttk::menubutton $f.first2c -menu $f.first2m  -width 16
        
        menu $f.first2m  -tearoff 0
	foreach i $mgsteps {
	    $f.first2m add command -label $mgsteplabel($i) -command "set meshoptions.firststep $i ; $f.first2c configure -text \"$mgsteplabel($i)\""
	}
        $f.first2m invoke $mgsteplabel(${meshoptions.firststep})        
        grid $f.first2l $f.first2c -sticky nw

        global meshoptions.laststep 
        ttk::label  $f.last2l -text "Last Step: "
        ttk::menubutton $f.last2c -menu $f.last2m -width 16

        menu $f.last2m  -tearoff 0

	foreach i $mgsteps {
	    $f.last2m add command -label $mgsteplabel($i) -command "set meshoptions.laststep $i ; $f.last2c configure -text \"$mgsteplabel($i)\""
	}
        $f.last2m invoke $mgsteplabel(${meshoptions.laststep})
        grid $f.last2l $f.last2c -sticky nw
        grid anchor $f center

	
	# tixOptionMenu $f.first -label "First Step : " \
	#     -options {
	# 	label.width  19
	# 	label.anchor e
	# 	menubutton.width 15
	#     } 

	# tixOptionMenu $f.last -label "Last Step : " \
	#     -options {
	# 	label.width  19
	# 	label.anchor e
	# 	menubutton.width 15
	#     } 

	# foreach step $mgsteps {
	#     $f.first add command $step -label $mgsteplabel($step)
	#     $f.last add command $step -label $mgsteplabel($step)
	# }

	# $f.first config -variable meshoptions.firststep 
	# $f.last config  -variable meshoptions.laststep 

	# pack $f.first $f.last
	




	set msg(0) "None"
	set msg(1) "Least"
	set msg(2) "Little"
	set msg(3) "Moderate"
	set msg(4) "Much"
	set msg(5) "Most"
	
	#tixOptionMenu $f.msg -label "Print Messages : " \
	    -options {
	#	label.width  19
	#	label.anchor e
	#	menubutton.width 15
	#    } 

        #foreach step {0 1 2 3 4 5 } {
        #    $f.msg add command $step -label $msg($step)
	#}
	#$f.msg config -variable options.printmsg 
	# pack $f.msg

        
        global options.printmsg
        #ttk::frame $f.msg2 
        ttk::label  $f.msg2l -text "Print Messages: "
        menu $f.msg2m  -tearoff 0
        ttk::menubutton $f.msg2c -menu $f.msg2m  -width 16
	foreach step {0 1 2 3 4 5 } {
	    $f.msg2m add command -label $msg($step) -command "set options.printmsg $step ; $f.msg2c configure -text $msg($step)"
            #            if { ${options.printmsg} == $step } { $f.msg2.c configure -text $msg($step) }
	}
        $f.msg2m invoke ${options.printmsg}
        grid $f.msg2l $f.msg2c -sticky nw
                
        set f $w.nb.general
        
        ttk::labelframe $f.bts -borderwidth 3 -relief groove -text "Additional meshing options"
        pack $f.bts -fill x -pady 15
	ttk::frame $f.bts.btnframe
	ttk::checkbutton $f.bts.btnframe.parthread -text "Parallel meshing thread" \
	    -variable options.parthread
	ttk::checkbutton $f.bts.btnframe.second -text "Second order elements" \
	    -variable options.secondorder
	ttk::checkbutton $f.bts.btnframe.quad -text "Quad dominated" \
	    -variable options.quad -command {
		if { ${options.quad} } {
		    set meshoptions.laststep os
		}
	    }
	ttk::checkbutton $f.bts.btnframe.invtets -text "Invert volume elements" \
	    -variable options.inverttets
	ttk::checkbutton $f.bts.btnframe.invtrigs -text "Invert surface elements" \
	    -variable options.inverttrigs
	ttk::checkbutton $f.bts.btnframe.azref -text "Automatic Z-refinement" \
	    -variable options.autozrefine
	pack $f.bts.btnframe -anchor center
	pack $f.bts.btnframe.parthread $f.bts.btnframe.second $f.bts.btnframe.quad $f.bts.btnframe.invtets $f.bts.btnframe.invtrigs $f.bts.btnframe.azref -anchor w


	# tixControl $f.elementorder -label "Element order: " -integer true \
	#     -variable options.elementorder -min 1 -max 20 \
	#     -options {
	# 	entry.width 2
	# 	label.width 20
	# 	label.anchor e
	#     }	
        # pack $f.elementorder

        #ttk::frame $f.bts.sbox        
        #pack $f.bts.sbox -anchor w -pady 10
        ttk::label $f.bts.btnframe.l -text "Element order"
        ttk::spinbox $f.bts.btnframe.elementorder2 -from 1 -to 20 -textvariable options.elementorder -width 2
        pack $f.bts.btnframe.elementorder2 $f.bts.btnframe.l  -anchor w -side left
        


        # ############################################################
        # Mesh - Size options
        # ############################################################        

        set f $w.nb.meshsize

        ttk::frame $f.f2
        pack $f.f2 -pady 10

        # # ttk::style configure Tframe -background red
        # puts "********************"
        # puts "found these themes:"
        # puts [ttk::themes]
        # ttk::setTheme classic
        # ttk::setTheme aqua
        # puts "style Tframe foreground = "
        # puts [ttk::style lookup Tframe -foreground]
        # puts "f2 style:"
        # puts [$f.f2 cget -style]
        # puts [winfo class $f.f2] 
        # puts "style element names gives:"
        # puts [ttk::style element names] 

        
        set f $f.f2
        
        ttk::frame $f.meshsize
        ttk::label $f.meshsize.l -text "max mesh-size"
        ttk::spinbox $f.meshsize.s -from 1e-9 -to 1e9 -textvariable options.meshsize -width 5 -validate focus -validatecommand "my_validatespinbox %W %P 10" \
	    -invalidcommand "my_invalidspinbox %W"
        pack $f.meshsize -fill x
        pack $f.meshsize.s $f.meshsize.l -side right

        ttk::frame $f.minmeshsize
        ttk::label $f.minmeshsize.l -text "min mesh-size"
        ttk::spinbox $f.minmeshsize.s -from 0 -to 1e9 -textvariable options.minmeshsize -width 5 -validate focus -validatecommand "my_validatespinbox %W %P 10" \
	    -invalidcommand "my_invalidspinbox %W"
        pack $f.minmeshsize -fill x
        pack $f.minmeshsize.s $f.minmeshsize.l -side right

        ttk::frame $f.grading
        ttk::label $f.grading.l -text "mesh-size grading"
        ttk::spinbox $f.grading.s -from 0.1 -to 1.0 -textvariable options.grading -width 5 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 3" \
            -invalidcommand "my_invalidspinbox %W"
        pack $f.grading -fill x
        pack $f.grading.s $f.grading.l -side right

        
	# tixControl $f.meshsize -label "max mesh-size: " -integer false \
	#     -variable options.meshsize -min 1e-9 -max 1e6 \
	#     -options {
	# 	entry.width 6
	# 	label.width 25
	# 	label.anchor e
	#     }	

	# tixControl $f.minmeshsize -label "min mesh-size: " -integer false \
	#     -variable options.minmeshsize -min 0 -max 1e6 \
	#     -options {
	# 	entry.width 6
	# 	label.width 25
	# 	label.anchor e
	#     }	

	# tixControl $f.grading -label "mesh-size grading: " -integer false \
	#     -variable options.grading -min 0.1 -max 1 -step 0.1 \
	#     -options {
	# 	entry.width 6
	# 	label.width 25
	# 	label.anchor e
	#     }	
	
	# pack $f.meshsize $f.minmeshsize $f.grading

        set f $w.nb.meshsize

	ttk::labelframe $f.msf -text "mesh-size file:" -relief groove -borderwidth 3
	pack $f.msf

 	# tixLabelEntry $f.msf.ent -label "mesh-size file: "  \
 	#     -labelside top \
 	#     -options {  
 	# 	entry.textVariable options.meshsizefilename 
 	# 	entry.width 35
 	# 	label.width 25
 	# 	label.anchor w
 	#     }

        ttk::entry $f.msf.ent -textvariable options.meshsizefilename -width 30
 	ttk::button $f.msf.btn -text "Browse" -command {
	    global options.meshsizefilename
	    set types {
		{"Meshsize file"   {.msz}	} }
	    set options.meshsizefilename [tk_getOpenFile -filetypes $types -initialfile ${options.meshsizefilename}]
	}

 	pack $f.msf.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $f.msf.btn -side left -anchor s -padx 4 -pady 4


	label $f.lab -text "Additional mesh size restrictions:"

	#csg-meshsize options

	ttk::labelframe $f.csg -relief groove -borderwidth 3 -text "CSG mesh-size"
	pack $f.csg -fill x
        proc test {a} {puts $a}
	#ttk::frame $f.csg.curv
	#pack $f.csg.curv -fill x -anchor center
        ttk::scale $f.csg.curvsc -orient horizontal -length 150 -from 0.2 -to 5 \
            -variable options.curvaturesafety -takefocus 0 -command "roundscale $f.csg.curvsc 1"
        #  -resolution 0.1 
        ttk::entry $f.csg.curve -textvariable options.curvaturesafety -width 3 \
            -validatecommand "my_validate %W [$f.csg.curvsc cget -from] [$f.csg.curvsc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
        ttk::label $f.csg.curvla -text "Elements per curvature radius"
	grid $f.csg.curvsc $f.csg.curve $f.csg.curvla -sticky nw -padx 4
        
	#ttk::frame $f.csg.elen 
	#pack $f.csg.elen  -fill x -anchor center
	ttk::scale $f.csg.elensc -orient horizontal -length 150 -from 0.2 -to 5 \
            -variable options.segmentsperedge -takefocus 0 -command "roundscale $f.csg.elensc 1"
        # -resolution 0.1
        ttk::entry $f.csg.elene -textvariable options.segmentsperedge -width 3 \
            -validatecommand "my_validate %W [$f.csg.elensc cget -from] [$f.csg.elensc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::label $f.csg.elenla -text "Elements per edge"
	grid $f.csg.elensc $f.csg.elene $f.csg.elenla -sticky nw -padx 4
        grid anchor $f.csg center
	

	#stl-meshsize options
	ttk::labelframe $f.stl -relief groove -borderwidth 3 -text "STL mesh-size"
	pack $f.stl -fill x 

	#ttk::frame $f.stl.r2
	#pack $f.stl.r2 -fill x
	ttk::scale $f.stl.r2sc -orient horizontal -length 150 -from 0.2 -to 5 \
            -variable stloptions.resthchartdistfac -takefocus 0 -command "roundscale $f.stl.r2sc 1"
        ttk::entry $f.stl.r2e -textvariable stloptions.resthchartdistfac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r2sc cget -from] [$f.stl.r2sc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::checkbutton $f.stl.r2bu -text "STL - chart distance" \
	    -variable stloptions.resthchartdistenable
	grid $f.stl.r2sc $f.stl.r2e $f.stl.r2bu -sticky nw -padx 4
	
	#ttk::frame $f.stl.r6
	#pack $f.stl.r6 -anchor w
	ttk::scale $f.stl.r6sc -orient horizontal -length 150 -from 0.2 -to 5 \
            -variable stloptions.resthlinelengthfac -takefocus 0 -command "roundscale $f.stl.r6sc 1"
        ttk::entry $f.stl.r6e -textvariable stloptions.resthlinelengthfac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r6sc cget -from] [$f.stl.r6sc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::checkbutton $f.stl.r6bu -text "STL - line length" \
	    -variable stloptions.resthlinelengthenable
	grid $f.stl.r6sc $f.stl.r6e $f.stl.r6bu -sticky nw -padx 4
	
	#ttk::frame $f.stl.r3
	#pack $f.stl.r3 -anchor w
	ttk::scale $f.stl.r3sc -orient horizontal -length 150 -from 0.2 -to 8 \
            -variable stloptions.resthcloseedgefac -takefocus 0 -command "roundscale $f.stl.r3sc 1"
        ttk::entry $f.stl.r3e -textvariable stloptions.resthcloseedgefac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r3sc cget -from] [$f.stl.r3sc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::checkbutton $f.stl.r3bu -text "STL/IGES/STEP - close edges" \
	    -variable stloptions.resthcloseedgeenable 
	
	grid $f.stl.r3sc $f.stl.r3e $f.stl.r3bu -sticky nw -padx 4
	
	#ttk::frame $f.stl.r1
	#pack $f.stl.r1 -anchor w
	ttk::scale $f.stl.r1sc -orient horizontal -length 150 -from 0.2 -to 5 \
	    -variable stloptions.resthsurfcurvfac -takefocus 0 -command "roundscale $f.stl.r1sc 1"
        ttk::entry $f.stl.r1e -textvariable stloptions.resthsurfcurvfac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r1sc cget -from] [$f.stl.r1sc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::checkbutton $f.stl.r1bu -text "STL - surface curvature" \
	    -variable stloptions.resthsurfcurvenable
	grid $f.stl.r1sc $f.stl.r1e $f.stl.r1bu -sticky nw -padx 4

	#ttk::frame $f.stl.r3b
	#pack $f.stl.r3b -anchor w
	ttk::scale $f.stl.r3bsc -orient horizontal -length 150 -from 0.2 -to 5 \
	    -variable stloptions.resthedgeanglefac -takefocus 0 -command "roundscale $f.stl.r3bsc 1"
        ttk::entry $f.stl.r3be -textvariable stloptions.resthedgeanglefac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r3bsc cget -from] [$f.stl.r3bsc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
        ttk::checkbutton $f.stl.r3bbu -text "STL - edge angle" \
            -variable stloptions.resthedgeangleenable
	grid $f.stl.r3bsc  $f.stl.r3be $f.stl.r3bbu -sticky nw -padx 4
	
	#ttk::frame $f.stl.r5
	#pack $f.stl.r5 -anchor w
	ttk::scale $f.stl.r5sc -orient horizontal -length 150 -from 0.2 -to 5 \
            -variable stloptions.resthsurfmeshcurvfac -takefocus 0 -command "roundscale $f.stl.r5sc 1"
        ttk::entry $f.stl.r5e -textvariable stloptions.resthsurfmeshcurvfac -width 3 \
            -validatecommand "my_validate %W [$f.stl.r5sc cget -from] [$f.stl.r5sc cget -to] %P 1" \
            -invalidcommand "my_invalid %W" -validate focus
	ttk::checkbutton $f.stl.r5bu -text "STL - surface mesh curv" \
	    -variable stloptions.resthsurfmeshcurvenable
	grid $f.stl.r5sc  $f.stl.r5e  $f.stl.r5bu -sticky nw -padx 4
	
	
	ttk::checkbutton $f.stl.recalch -text "STL - Recalc mesh size for surface optimization" \
	    -variable stloptions.recalchopt
	grid $f.stl.recalch -sticky n -columnspan 3 -column 0

	ttk::button $f.stl.calch -text "Calc New H" -command { redraw; Ng_STLCalcLocalH }
	grid $f.stl.calch -columnspan 3 -column 0
	grid anchor $f.stl center
        # set f [$w.nb subwidget chartopt]
	
        # round ttk::scale values to n_digits
        proc roundscale {w n_digits args} {
            set val [$w get]
	    global [$w cget -variable]
	    if {$n_digits == 0 } {
                set [$w cget -variable] [tcl::mathfunc::round $val]
	    } else {
                set [$w cget -variable] [format "%.[append n_digits "f"]" $val]
	    }
	}

        # validate ttk::entry which are linked to ttk::scales widgets
        global last_accepted_sc
        proc my_validate {w mini maxi val n_digits} {
            global last_accepted_sc [$w cget -textvariable]            
            if {[string length $val] == 0} {return 0}
            if {[string is double $val] == 1} {
                if { $n_digits == 0 } {
                    set val [tcl::mathfunc::max $mini [tcl::mathfunc::min $maxi [tcl::mathfunc::round $val]]]                    
                } else {
                    if { $n_digits < 9 } {
                        set val [tcl::mathfunc::max $mini [tcl::mathfunc::min $maxi [format "%.[append n_digits "f"]" $val]]]
                    }
                }                    
                    set last_accepted_sc $val
                    set [$w cget -textvariable] $val
                    return 1
                } else {
                    return 0
                }
        }

        # if my_validate returns 0, this function gets called
        proc my_invalid {w} {
            global last_accepted_sc [$w cget -textvariable]
            set [$w cget -textvariable] $last_accepted_sc
        }

        set f $w.nb.chartopt
        ttk::labelframe $f.mainframe -text "STL angles" -relief groove -borderwidth 3

        pack $f.mainframe -fill x -pady 15
        set f $f.mainframe

        #ttk::frame $f.f1
        ttk::label $f.labYangles -text "Yellow Edges Angle ()"
        ttk::scale $f.scale1 -orient horizontal -length 150 -from 0 -to 90 -variable stloptions.yangle -takefocus 0 -command "roundscale $f.scale1 1"
        ttk::entry $f.entry1 -textvariable stloptions.yangle -width 5 -validate focus -takefocus 0 -validatecommand "my_validate %W [$f.scale1 cget -from] [$f.scale1 cget -to] %P 1" \
            -invalidcommand "my_invalid %W"

        #pack $f.f1 -anchor center
        grid $f.scale1 $f.entry1 $f.labYangles -sticky nw -padx 4 -pady 6

        #ttk::frame $f.f21
        ttk::label $f.labEangles -text "Edge Corner Angle ()"
        ttk::scale $f.scale2 -orient horizontal -length 150 -from 0 -to 180 -variable stloptions.edgecornerangle -takefocus 0 -command "roundscale $f.scale2 1"
        ttk::entry $f.entry2 -textvariable stloptions.edgecornerangle -width 5 -validate focus -takefocus 0 -validatecommand "my_validate %W [$f.scale2 cget -from] [$f.scale2 cget -to] %P 1" \
            -invalidcommand "my_invalid %W"

        #pack $f.f21 -anchor center
        grid $f.scale2 $f.entry2 $f.labEangles -sticky nw -padx 4 -pady 6

        #ttk::frame $f.f31
        ttk::label $f.lab31 -text "Chart Angle ()"
        ttk::scale $f.scale3 -orient horizontal -length 150 -from 0 -to 180 -variable stloptions.chartangle -takefocus 0 -command "roundscale $f.scale3 1"
        ttk::entry $f.entry3 -textvariable stloptions.chartangle -width 5 -validate focus -takefocus 0 -validatecommand "my_validate %W [$f.scale3 cget -from] [$f.scale3 cget -to] %P 1" \
            -invalidcommand "my_invalid %W"

        #pack $f.f31 -anchor center
        grid $f.scale3 $f.entry3 $f.lab31 -sticky nw -padx 4 -pady 6

        #ttk::frame $f.f41
        ttk::label $f.lab41 -text "Outer Chart Angle ()"
        ttk::scale $f.scale4 -orient horizontal -length 150 -from 0 -to 180 -variable stloptions.outerchartangle -takefocus 0 -command "roundscale $f.scale4 1"
        ttk::entry $f.entry4 -textvariable stloptions.outerchartangle -width 5 -validate focus -takefocus 0 -validatecommand "my_validate %W [$f.scale4 cget -from] [$f.scale4 cget -to] %P 1" \
            -invalidcommand "my_invalid %W"

        #pack $f.f41 -anchor center
        grid $f.scale4 $f.entry4 $f.lab41 -sticky nw -padx 4 -pady 6
	grid anchor $f center
	# Optimization options

        global last_accepted_sp
        # Used to validate the entries linked with a ttk::spinbox widget
        proc my_validatespinbox {w val n_digits} {
            global last_accepted_sp            
                if {[string length $val] == 0} {return 0}
                if {[string is double $val] == 1} {
                    if { $n_digits == 0 } {
                        if { $n_digits < 9 } {
                            set val [tcl::mathfunc::round $val] } else { set val [format "%.[append n_digits "f"]" $val]
                        }
                    }
                        $w set [tcl::mathfunc::max [$w cget -from] [tcl::mathfunc::min [$w cget -to] $val]]
                        set last_accepted_sp $val                        
                        return 1
                } else {
                    return 0
                }
        }
	
        proc my_invalidspinbox {w} {
            global last_accepted_sp
            $w set $last_accepted_sp
        }
        # set f [$w.nb subwidget optimizer]
        set f $w.nb.optimizer
        ttk::labelframe $f.optframe -text "Optimization settings" -relief groove -borderwidth 3
	pack $f.optframe -fill x -pady 15
	
        #ttk::frame $f.optframe.sos
        ttk::label $f.optframe.sosl -text "Surface opt steps"
        ttk::spinbox $f.optframe.soss -from 0 -to 99 -textvariable options.optsteps2d -width 5 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        #pack $f.optframe.sos -anchor center
        grid $f.optframe.sosl $f.optframe.soss -sticky nw;# -side right -fill x -pady 2

        #ttk::frame $f.optframe.vos
        ttk::label $f.optframe.vosl -text "Volume opt steps"
        ttk::spinbox $f.optframe.voss -from 0 -to 99 -textvariable options.optsteps3d -width 5 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        #pack $f.optframe.vos -anchor center
        grid $f.optframe.vosl $f.optframe.voss -sticky nw;# -side right -fill x -pady 2
	
        #ttk::frame $f.optframe.esw
        ttk::label $f.optframe.eswl -text "Element size weight"
        ttk::spinbox $f.optframe.esws -from 0 -to 1 -textvariable options.elsizeweight -width 5 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 1" \
            -invalidcommand "my_invalidspinbox %W"

        #pack $f.optframe.esw -anchor center
        grid $f.optframe.eswl $f.optframe.esws -sticky nw;# -side right -fill x -pady 2

        #ttk::frame $f.optframe.wem
        ttk::label $f.optframe.weml -text "Worst element measure"
        ttk::spinbox $f.optframe.wems -from 1 -to 10 -textvariable options.opterrpow -width 5 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        #pack $f.optframe.wem -anchor e
        grid $f.optframe.weml $f.optframe.wems -sticky nw;# -side right -fill x -pady 2
        grid anchor $f.optframe center 
        # These functions are needed due to a bug within the aqua theme
        # if a ttk::scale widget has a from value larger than 100.
	proc roundscale_helper_osx {w val} {
            global [$w cget -variable] options.badellimit
            set [$w cget -variable] [tcl::mathfunc::round $val]
            set options.badellimit [expr [tcl::mathfunc::round $val]+160]
        }

        proc my_validate_helper_osx {w val} {
            if {[string length $val] == 0} {return 0}
            if {[string is double $val] == 1} {			
                set scale_loc [lindex [winfo children [winfo parent $w]] [lsearch [winfo children [winfo parent $w]] *scale]]            
                global [$scale_loc cget -variable] options.badellimit
                set [$scale_loc cget -variable] [tcl::mathfunc::max [$scale_loc cget -from] [tcl::mathfunc::min [$scale_loc cget -to] [expr [tcl::mathfunc::round $val]-160]]]
                set options.badellimit [tcl::mathfunc::max [expr [$scale_loc cget -from]+160] [tcl::mathfunc::min [expr [$scale_loc cget -to]+160] [tcl::mathfunc::round $val]]]
                return 1
            } else {
                return 0
            }
        }

        proc my_invalid_helper_osx {w} {
            global options.badellimit
            set scale_loc [lindex [winfo children [winfo parent $w]] [lsearch [winfo children [winfo parent $w]] *scale]]
            global [$scale_loc cget -variable]
            set [$scale_loc cget -variable] [tcl::mathfunc::round [$scale_loc get]]
            set options.badellimit [expr [tcl::mathfunc::round [$scale_loc get]]+160]
        }    
    
        global dummy_badellimit
        set dummy_badellimit 15 
        ttk::labelframe $f.optframe2 -text "Bad elements" -relief groove -borderwidth 3
        pack $f.optframe2 -fill x -pady 15 -ipady 5
        ttk::frame $f.optframe2.badellimit        
        ttk::label $f.optframe2.lab -text "bad element criterion";
        ttk::scale $f.optframe2.scale -orient horizontal -length 100 -from 00 -to 20 -variable dummy_badellimit -takefocus 0 -command "roundscale_helper_osx $f.optframe2.scale"
        ttk::entry $f.optframe2.entry -textvariable options.badellimit -width 3 -validate focusout -takefocus 0 -validatecommand "my_validate_helper_osx %W %P" \
            -invalidcommand "my_invalid_helper_osx %W"
        #pack $f.optframe2.badellimit -anchor center
        grid $f.optframe2.scale $f.optframe2.entry $f.optframe2.lab -padx 4 -sticky nw
        grid anchor $f.optframe2 center



	# insider options
        # set f [$w.nb subwidget insider]
        set f $w.nb.debug    
        ttk::labelframe $f.f2 -text "Advanced options" -borderwidth 3 -relief groove
        pack $f.f2 -fill x -pady 15
        #ttk::frame $f.f2.frame
        #pack $f.f2.frame 
        set f $f.f2
	ttk::checkbutton $f.localh -text "Use Local Meshsize" \
	    -variable options.localh
	ttk::checkbutton $f.delauney -text "Use Delaunay" \
	    -variable options.delaunay
	ttk::checkbutton $f.checkoverlap -text "Check Overlapping" \
	    -variable options.checkoverlap
	ttk::checkbutton $f.checkcb -text "Check Chart Boundary" \
	    -variable options.checkchartboundary
	ttk::checkbutton $f.blockfill -text "Do Blockfilling" \
	    -variable options.blockfill

        grid $f.localh  $f.delauney -sticky nw
        grid $f.checkoverlap $f.blockfill -sticky nw
        grid $f.checkcb -sticky nw
	grid anchor $f center
        # debugging options    
        set f $w.nb.debug

        # enable / disable ttk::entry widgets linked to ttk::checkbuttons
        proc enable_cb {w1 w2 w3} {
            Ng_SetDebugParameters
            if {[string match *selected* [$w1 state]] == 1 } {
                $w2 configure -state normal
                $w3 configure -state normal
            } else {
                $w2 configure -state disabled
                $w3 configure -state disabled
            }
        }
	
        ttk::labelframe $f.cb1 -text "Debugging options" -borderwidth 3 -relief groove
        pack $f.cb1 -fill x -pady 15
        
        #frame $f.cb1.cb0
        #pack $f.cb1.cb0 -fill x 
        
        ttk::checkbutton $f.cb1.slowchecks -text "Slow checks" \
            -variable debug.slowchecks -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.debugoutput -text "Debugging outout" \
            -variable debug.debugoutput -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltexline -text "Halt on exising line" \
            -variable debug.haltexistingline  -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltoverlap -text "Halt on Overlap" \
            -variable debug.haltoverlap  -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltsuc -text "Halt on success" \
            -variable debug.haltsuccess  -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltnosuc -text "Halt on no success" \
            -variable debug.haltnosuccess  -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltlargequal -text "Halt on large quality class" \
            -variable debug.haltlargequalclass  -command { Ng_SetDebugParameters }
        ttk::checkbutton $f.cb1.haltseg -text "Halt on Segment:" \
            -variable debug.haltsegment  -command "enable_cb %W $f.cb1.segs.ent1 $f.cb1.segs.ent2"
        ttk::checkbutton $f.cb1.haltnode -text "Halt on Node:" \
            -variable debug.haltnode  -command "enable_cb %W $f.cb1.segs.ent1 $f.cb1.segs.ent2"
        ttk::frame $f.cb1.fr
        ttk::checkbutton $f.cb1.fr.cb -text "Halt on Face:" \
            -variable debug.haltface -command "enable_cb %W $f.cb1.fr.ent $f.cb1.fr.ent"
        ttk::entry $f.cb1.fr.ent -textvariable debug.haltfacenr -width 3 -state disabled
        
        pack $f.cb1.fr.cb $f.cb1.fr.ent -side left 

        ttk::frame $f.cb1.segs
        ttk::label $f.cb1.segs.lab1 -text "P1:"
        ttk::entry $f.cb1.segs.ent1 -width 6 \
            -textvariable debug.haltsegmentp1  -state disabled
        ttk::label $f.cb1.segs.lab2 -text "P2:"
        ttk::entry $f.cb1.segs.ent2 -width 6 \
            -textvariable debug.haltsegmentp2  -state disabled

        pack $f.cb1.segs.lab1 $f.cb1.segs.ent1 $f.cb1.segs.lab2 $f.cb1.segs.ent2 -side left


        grid $f.cb1.slowchecks $f.cb1.debugoutput -sticky nw
        grid $f.cb1.haltexline $f.cb1.haltoverlap -sticky nw
        grid $f.cb1.haltsuc $f.cb1.haltnosuc  -sticky nw
        grid $f.cb1.haltlargequal  $f.cb1.fr -sticky nw
        grid $f.cb1.haltnode -sticky nw
        grid $f.cb1.haltseg  -stick nw
        grid $f.cb1.segs -stick w -row 4 -rowspan 2 -column 1

        grid rowconfigure $f.cb1 3 -pad 8
        grid anchor $f.cb1 center
        ttk::checkbutton $f.cb1.showactivechart -text "Show Active Meshing-Chart" -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }

        grid $f.cb1.showactivechart
        grid rowconfigure $f.cb1 3 -pad 8
        grid rowconfigure $f.cb1 5 -pad 8
        
        set f $w.nb.debug        
        ttk::labelframe $f.cont -relief groove -borderwidth 3 -text "Debugging visualization"
        pack $f.cont -fill x -pady 15
        #ttk::frame $f.cont.f
        #pack $f.cont.f

        ttk::checkbutton $f.cont.multidrawing -text "Draw Meshing" -variable multithread_drawing
        ttk::checkbutton $f.cont.multitestmode -text "Meshing Testmode" -variable multithread_testmode
        ttk::button $f.cont.goon -text "Go On" -command { set multithread_pause 0 }

        grid $f.cont.multidrawing -sticky nw
        grid $f.cont.multitestmode -sticky nw
        grid $f.cont.goon -row 0 -rowspan 2 -column 1 -sticky w
        grid columnconfigure $f.cont 0 -pad 30
        grid columnconfigure $f.cont 1 -pad 20
        grid anchor $f.cont center

	global userlevel
	if { $userlevel < 3} {
	    $w.nb delete insider
	    $w.nb delete debug
	}



# 	tixButtonBox $w.bbox -orientation horizontal
# 	$w.bbox add ok    -text Apply  -underline 0 -width 5 \
# 	    -command { 
# 		[.options_dlg.nb subwidget meshsize].meshsize invoke
# 		[.options_dlg.nb subwidget meshsize].grading invoke
# 		[.options_dlg.nb subwidget optimizer].os2d invoke
# 		[.options_dlg.nb subwidget optimizer].os3d invoke
# 		[.options_dlg.nb subwidget optimizer].elw invoke
# 		[.options_dlg.nb subwidget optimizer].wem invoke
		
# 		Ng_SetMeshingParameters 
# 	    }
	
# 	$w.bbox add close -text Done   -underline 0 -width 5 \
# 	    -command {
# 		[.options_dlg.nb subwidget meshsize].meshsize invoke
# 		[.options_dlg.nb subwidget meshsize].grading invoke
# 		[.options_dlg.nb subwidget optimizer].os2d invoke
# 		[.options_dlg.nb subwidget optimizer].os3d invoke
# 		[.options_dlg.nb subwidget optimizer].elw invoke
# 		[.options_dlg.nb subwidget optimizer].wem invoke
		
# 		Ng_SetMeshingParameters
# 		destroy .options_dlg
# 	    }
	
#       pack $w.bbox -side bottom -fill x


        ttk::frame $w.bu
        pack $w.bu -fill x -ipady 3

        ttk::button $w.bu.apl -text "Apply" -command { 	    
            Ng_SetMeshingParameters 
            Ng_SetDebugParameters
        }

        ttk::button $w.bu.ok -text "Done" -command {
            Ng_SetMeshingParameters
            Ng_SetDebugParameters
            wm withdraw .options_dlg
#           destroy .options_dlg
        }

        pack  $w.bu.apl $w.bu.ok -side left -expand yes 
        wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Meshing Options"
	focus .options_dlg
    }
}


meshingoptionsdialog
wm withdraw .options_dlg









#
#
#  Viewing dialog
#
#
proc viewingoptionsdialog { } {

    global userlevel

    set w .viewopts_dlg
    
    if {[winfo exists .viewopts_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w
        #wm resizable $w 0 0 
 
        
        pack [ttk::notebook $w.nb]  -fill both -fill both -side top
        $w.nb add [ttk::frame $w.nb.general] -text "General" -underline 0
        $w.nb add [ttk::frame $w.nb.stl] -text "STL" -underline 0
        $w.nb add [ttk::frame $w.nb.occ] -text "IGES/STEP" -underline 0
        $w.nb add [ttk::frame $w.nb.mesh] -text "Mesh" -underline 0
        $w.nb add [ttk::frame $w.nb.light] -text "Light" -underline 0
        $w.nb add [ttk::frame $w.nb.edges] -text "Edges" -underline 0
        $w.nb add [ttk::frame $w.nb.misc] -text "Misc." -underline 3
        

	# tixNoteBook $w.nb -ipadx 6 -ipady 6
	
	# $w.nb add general -label "General" -underline 0
	# $w.nb add stl -label "STL" -underline 0
	# $w.nb add occ -label "IGES/STEP" -underline 0
	# $w.nb add mesh -label "Mesh"   -underline 0
	# $w.nb add light -label "Light"   -underline 0
	# $w.nb add edges -label "Edges"   -underline 0
	# $w.nb add misc -label "Misc."  -underline 3

	# pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	





	# general
	set f $w.nb.general
        ttk::labelframe $f.gvop -text "General viewing options" -relief groove -borderwidth 3
        pack $f.gvop -fill x -pady 15
        set f $f.gvop
	ttk::checkbutton $f.backcol -text "White Background" \
	-variable viewoptions.whitebackground \
	-command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.cross -text "Draw Coordinate Cross" \
	-variable viewoptions.drawcoordinatecross \
	-command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.color -text "Draw Color-bar" \
	-variable viewoptions.drawcolorbar \
	-command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.netgen -text "Draw Netgen-logo" \
	-variable viewoptions.drawnetgenlogo \
	-command { Ng_SetVisParameters; redraw }

	grid $f.backcol -sticky nw
        grid $f.cross -stick nw
        grid $f.color  -sticky nw
        grid $f.netgen -sticky nw
# 	checkbutton $f.stereo -text "Stereo View" \
# 	-variable viewoptions.stereo \
# 	-command { Ng_SetVisParameters; redraw }
# 	pack $f.stereo


        menu $f.stylemenu
        ttk::menubutton $f.style -menu $f.stylemenu -width 10 -text [ttk::style theme use]
        grid $f.style -sticky nw 
        grid anchor $f center
        
        foreach theme [ttk::themes] {
            $f.stylemenu add command -label  $theme \
                -command " $f.style configure -text $theme; puts $theme ; ttk::setTheme $theme"
        }
        
	# stl geometry 
	set f $w.nb.stl
	ttk::labelframe $f.show -relief groove -borderwidth 3 -text "STL viewing options"
	pack $f.show -fill x -pady 15
	ttk::checkbutton $f.show.showtrias -text "Show STL-Triangles" \
	    -variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }
	#grid $f.show.showtrias -stick nw
	
	ttk::checkbutton $f.show.showfilledtrias -text "Show Filled Triangles" \
	    -variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }
	grid $f.show.showtrias $f.show.showfilledtrias -sticky nw
	
	ttk::checkbutton $f.show.showactivechart -text "Show Active Meshing-Chart" \
	    -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }
	#grid $f.show.showactivechart -sticky nw
	
	ttk::checkbutton $f.show.showedges -text "Show Edges" \
	    -variable stloptions.showedges -command { Ng_SetVisParameters; redraw }
	grid $f.show.showactivechart $f.show.showedges -sticky nw
	grid anchor $f.show center
	#frame $f.special -relief groove -borderwidth 3
	#pack $f.special
	ttk::checkbutton $f.show.showmarktrias -text "Show Chart Triangles" \
	    -variable stloptions.showmarktrias \
	    -command {set stldoctor.showfaces 0; Ng_STLDoctor; Ng_SetVisParameters; redraw }
	#pack $f.show.showmarktrias -side left

	ttk::checkbutton $f.show.showfaces -text "Show Faces" \
	    -variable stldoctor.showfaces \
	    -command {set stloptions.showmarktrias 0; Ng_STLDoctor; Ng_SetVisParameters; redraw}    
	#pack $f.show.showfaces -side left
        grid $f.show.showmarktrias $f.show.showfaces -sticky nw

        
        ttk::labelframe $f.fn -relief groove -borderwidth 3 -text "Chart/Face number"
	pack $f.fn -fill x 
	ttk::label $f.fn.lab3 -text "Chart/Face number"        
	ttk::scale $f.fn.scale3 -orient horizontal -length 150 -from 0 -to 200 \
            -variable stloptions.chartnumber -command "Ng_SetVisParameters; redraw;roundscale $f.fn.scale3 0"
        ttk::entry $f.fn.ent3 -textvariable stloptions.chartnumber -width 3 -validate focus -takefocus 0 \
            -validatecommand "Ng_SetVisParameters; redraw;my_validate %W [$f.fn.scale3 cget -from] [$f.fn.scale3 cget -to] %P 0" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"
            
	grid $f.fn.scale3 $f.fn.ent3 $f.fn.lab3 -sticky nw -padx 4        
        
	
	#frame $f.fo -relief groove -borderwidth 3
	#pack $f.fo
	tk::label $f.fn.lab -text "Chart/Face Offset:";
	ttk::entry $f.fn.ent -width 3 \
	    -textvariable stloptions.chartnumberoffset -validate focus -takefocus 0 \
            -validatecommand "my_validate %W 0 1e9 %P 0" \
            -invalidcommand "my_invalid %W"
	grid $f.fn.lab -sticky ne -padx 4
        grid $f.fn.ent -sticky nw -padx 4 -row 1 -column 1
        grid anchor $f.fn center 
        
	ttk::labelframe $f.advstl -text "Advanced STL options" -relief groove -borderwidth 3
        pack $f.advstl -fill x -pady 15
        #frame $f.mt
	#pack $f.mt -fill x
	ttk::checkbutton $f.advstl.bu1 -text "Show Marked (Dirty) Triangles" \
	    -variable stldoctor.showmarkedtrigs \
	    -command {Ng_STLDoctor; redraw}    
	#pack $f.mt.bu

	#frame $f.ep
	#pack $f.ep -fill x
	ttk::checkbutton $f.advstl.bu2 -text "show edge corner points" \
	    -variable stldoctor.showedgecornerpoints \
	    -command {Ng_STLDoctor; redraw}    
	#pack $f.ep.bu

	#frame $f.stt
	#pack $f.stt -fill x
	ttk::checkbutton $f.advstl.bu3 -text "show touched triangle chart" \
	    -variable stldoctor.showtouchedtrigchart \
	    -command {set stldoctor.showfaces 0; set stloptions.showmarktrias 1; \
			  Ng_STLDoctor; Ng_SetVisParameters; redraw}    
	#pack $f.stt.bu

	#frame $f.sml
	#pack $f.sml -fill x
	ttk::checkbutton $f.advstl.bu4 -text "draw meshed edges" \
	    -variable stldoctor.drawmeshededges \
	    -command {Ng_STLDoctor;}    
	#pack $f.sml.bu
	
	
	#frame $f.sm
	#pack $f.sm -fill x
	ttk::checkbutton $f.advstl.bu5 -text "select with mouse" \
	    -variable stldoctor.selectwithmouse
	#pack $f.sm.bu
	grid $f.advstl.bu1 -stick nw
        grid $f.advstl.bu2 -sticky nw
        grid $f.advstl.bu3 -stick nw
        grid $f.advstl.bu4 -stick nw
        grid $f.advstl.bu5 -stick nw
        grid anchor $f.advstl center
	ttk::frame $f.advstl.tbn
	ttk::label $f.advstl.tbn.lab -text "Select triangle by number";
	ttk::entry $f.advstl.tbn.ent -width 5 \
	    -textvariable stldoctor.selecttrig 
        pack $f.advstl.tbn.lab $f.advstl.tbn.ent -padx 4 -side left
	grid $f.advstl.tbn -sticky nw
	grid anchor $f.advstl center
        grid rowconfigure $f.advstl 4 -pad 8
	
        ttk::labelframe $f.vc -relief groove -borderwidth 3 -text "Vicinity options"
	pack $f.vc -fill x -pady 15
	ttk::checkbutton $f.vc.bu -text "show vicinity" \
	    -variable stldoctor.showvicinity \
	    -command {Ng_STLDoctor vicinity; redraw}
	ttk::label $f.vc.lab -text "vicinity size";
	ttk::scale $f.vc.scale -orient horizontal -length 150 -from 0 -to 200 \
	    -variable stldoctor.vicinity \
             -takefocus 0 \
             -command "roundscale $f.vc.scale 0; Ng_STLDoctor vicinity; redraw"
	    #-command { Ng_STLDoctor vicinity; redraw }
        ttk::entry $f.vc.ent -width 4 -textvariable stldoctor.vicinity -validate focus \
        -takefocus 0 -validatecommand "Ng_STLDoctor vicinity; redraw;my_validate %W [$f.vc.scale cget -from] [$f.vc.scale cget -to] %P 0" \
            -invalidcommand "my_invalid %W;Ng_STLDoctor vicinity; redraw"
	#pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes
	grid $f.vc.bu -stick nw -columnspan 3 -column 0 
        grid $f.vc.scale $f.vc.ent $f.vc.lab -sticky nw -padx 4
        grid anchor $f.vc center


	# IGES/STEP
	set f $w.nb.occ
	ttk::labelframe $f.occframe -text "IGES/STEP options" -relief groove -borderwidth 3
        pack $f.occframe -fill x -pady 15 -ipady 8
        #set f $f.occframe
	ttk::checkbutton $f.occframe.occshowsurfaces -text "Show surfaces " \
	    -variable occoptions.showsurfaces \
	    -command { Ng_SetOCCVisParameters; redraw }

	ttk::checkbutton $f.occframe.occshowedges -text "Show edges " \
	    -variable occoptions.showedges \
	    -command { Ng_SetOCCVisParameters; redraw }
        grid $f.occframe.occshowsurfaces $f.occframe.occshowedges -sticky nw -padx 4
        grid anchor $f.occframe center
	#ttk::frame $f.deflection -relief groove -borderwidth 3
	#pack $f.occframe.deflection -fill x
	ttk::button $f.occframe.btn -text "Rebuild visualization data" \
	    -command {
		Ng_SetOCCVisParameters
		Ng_OCCCommand buildvisualizationmesh
		redraw
	    }

	#tixControl $f.occframe.ent -label "Visualization smoothness" -integer false \
	    -variable occoptions.deflection -min 0.1 -max 3 -step 0.1 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters }
        ttk::frame $f.occframe.vssm
        ttk::label $f.occframe.vssm.lab -text "Visulization smoothness"        
        ttk::spinbox $f.occframe.vssm.sp -textvariable occoptions.deflection \
            -from 0.1 -to 3 -increment 0.1 -width 4 -command { catch Ng_SetOCCVisParameters } \
            -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"
        pack $f.occframe.vssm.lab $f.occframe.vssm.sp -side left -padx 4
        grid $f.occframe.vssm -sticky nw -columnspan 2 -column 0 -pady 8        
        grid $f.occframe.btn -columnspan 2 -column 0 -sticky n



	#grid $f.occframe.ent $f.occframe.lab -sticky nw	


	# ACIS visualization / construction
        
        ttk::labelframe $f.occframe1 -relief groove -borderwidth 3 -text "ACIS visulization / construction"
        pack $f.occframe1 -fill x -pady 15 -ipady 8 
        #ttk::frame $f.occframe1.shso
        ttk::label $f.occframe1.lab1 -text "Show solid (0 for all)"
        ttk::spinbox $f.occframe1.sp1 -textvariable occoptions.showsolidnr \
            -from 0 -to 999 -increment 1 -width 4 -command { catch Ng_SetOCCVisParameters;redraw } \
            -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"
        #pack $f.occframe1.shso.lab $f.occframe1.shso.sp -side left -padx 4

        #ttk::frame $f.occframe1.shso2
        ttk::label $f.occframe1.lab2 -text "Show solid 2"
        ttk::spinbox $f.occframe1.sp2 -textvariable occoptions.showsolidnr2 \
            -from 0 -to 999 -increment 1 -width 4 -command { catch Ng_SetOCCVisParameters;redraw } \
            -validate focus -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"
        #pack $f.occframe1.shso2.lab $f.occframe1.shso2.sp -side left -padx 4
        
	#tixControl $f.showsolid -label "Show solid (0 for all)" -integer true \
            -variable occoptions.showsolidnr -min 0 -max 999 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters; redraw }
    
	#tixControl $f.showsolid2 -label "Show solid 2" -integer true \
            -variable occoptions.showsolidnr2 -min 0 -max 999 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters; redraw }

	ttk::button $f.occframe1.subtract -text "Subtract (2 minus 1)" \
	    -command {
		Ng_ACISCommand subtract ${occoptions.showsolidnr} ${occoptions.showsolidnr2}
		redraw
	    }

            
            
	ttk::button $f.occframe1.combine -text "Combine all" \
	    -command {
		Ng_ACISCommand combineall
		redraw
	    }

	#pack $f.showsolid $f.showsolid2 $f.subtract $f.combine;# -sticky nw
        grid $f.occframe1.lab1 -row 0 -column 0 -sticky ne
        grid $f.occframe1.sp1 -row 0 -column 1 -sticky nw
        grid $f.occframe1.lab2 -row 1 -column 0 -sticky ne
        grid $f.occframe1.sp2 -row 1 -column 1 -sticky nw
        grid $f.occframe1.combine -columnspan 2 -column 0 -sticky n
        grid anchor $f.occframe1 center        




	# mesh options
	set f $w.nb.mesh
        
	ttk::labelframe $f.center -relief groove -borderwidth 3 -text "how shall i name you?"
	pack $f.center -fill x -pady 15
	ttk::button $f.center.lab1 -text "Set Center Point" \
	    -command { Ng_SetVisParameters; Ng_Center; redraw }
	ttk::entry $f.center.ent1 -width 5 \
	    -textvariable viewoptions.centerpoint -validate focus \
            -validatecommand "my_validate %W 0 1e9 %P 0" \
            -invalidcommand "my_invalid %W"
	grid $f.center.ent1 $f.center.lab1 -padx 4 -pady 4 -sticky nw
	
	#ttk::frame $f.drawel -relief groove -borderwidth 3
	#pack $f.drawel -fill x
	ttk::button $f.center.lab2 -text "Draw Element" \
	    -command { Ng_SetVisParameters; Ng_ZoomAll; redraw }
	ttk::entry $f.center.ent2 -width 5 \
	    -textvariable viewoptions.drawelement -validate focus \
            -validatecommand "my_validate %W 0 1e9 %P 0" \
            -invalidcommand "my_invalid %W"
	grid $f.center.ent2 $f.center.lab2 -padx 4 -pady 4 -sticky nw
        grid anchor $f.center center
        
        ttk::labelframe $f.meshframe -text "Mesh visualization options" -relief groove -borderwidth 3
        pack $f.meshframe -fill x -pady 15
        set f $f.meshframe
	ttk::checkbutton $f.showcolor -text "Meshsize Visualization" \
	    -variable viewoptions.colormeshsize \
	    -command { Ng_SetVisParameters;redraw; }

	ttk::checkbutton $f.showfilledtrigs -text "Show filled triangles" \
	    -variable viewoptions.drawfilledtrigs \
	    -command { Ng_SetVisParameters; redraw }
	
	ttk::checkbutton $f.showedges -text "Show edges" \
	    -variable viewoptions.drawedges \
	    -command { Ng_SetVisParameters; redraw }
	
	ttk::checkbutton $f.showoutline -text "Show Triangle Outline" \
	    -variable viewoptions.drawoutline \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showbadels -text "Show bad elements" \
	    -variable viewoptions.drawbadels \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showprisms -text "Show prisms" \
	    -variable viewoptions.drawprisms \
	    -command { Ng_SetVisParameters; redraw }
            
	ttk::checkbutton $f.showpyramids -text "Show pyramids" \
	    -variable viewoptions.drawpyramids \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showhexes -text "Show hexes" \
	    -variable viewoptions.drawhexes \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showidentified -text "Show identified points" \
	    -variable viewoptions.drawidentified \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showmetispartition -text "Show METIS Partition" \
	    -variable viewoptions.drawmetispartition \
	    -command { Ng_SetVisParameters; redraw }

	ttk::checkbutton $f.showpointnumbers -text "Show Point-numbers" \
	    -variable viewoptions.drawpointnumbers \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showedgenumbers -text "Show Edge-numbers" \
	    -variable viewoptions.drawedgenumbers \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showfacenumbers -text "Show Face-numbers" \
	    -variable viewoptions.drawfacenumbers \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showelementnumbers -text "Show Element-numbers" \
	    -variable viewoptions.drawelementnumbers \
	    -command { Ng_SetVisParameters; redraw }
	
	# label $f.showdomainlab -text "Domain Surface"
#	scale $f.showdomain -orient horizontal -length 100 -from 0 -to 50 \
	    -resolution 1 -variable  viewoptions.drawdomainsurf    \
	    -command { Ng_SetVisParameters; redraw } \
	    -label "Domain Surface" 

        #pack $f.showfilledtrigs
	#pack $f.showoutline $f.subdiv $f.showedges  $f.showbadels 
	## pack $f.showdomainlab 
	#pack $f.showdomain 
	#pack $f.showpointnumbers 
	#pack $f.showedgenumbers $f.showfacenumbers $f.showelementnumbers 
	#pack $f.showmetispartition


        
        
	ttk::frame $f.frametets
	ttk::checkbutton $f.frametets.showtets -text "" \
	    -variable viewoptions.drawtets \
	    -command { Ng_SetVisParameters; redraw }
        ttk::label $f.frametets.label -text "\Show Tets\rin domain"
        ttk::spinbox $f.frametets.showtetsdomain -from 0 -to 500 -increment 1 -width 3 \
            -textvariable viewoptions.drawtetsdomain -validate focus \
            -command "Ng_SetVisParameters; redraw;" \
            -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"
        
        #ttk::frame $f.frametets
        ttk::label $f.frametets.label1 -text "Subdivision"
        ttk::spinbox $f.frametets.subdiv -from 0 -to 8 -increment 1 -width 3 \
            -textvariable visoptions.subdivisions -validate focus \
            -command { Ng_SetVisParameters; Ng_Vis_Set parameters; Ng_SetNextTimeStamp; redraw } \
            -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        ttk::label $f.frametets.label2 -text "Show surface\rof domain"
        ttk::spinbox $f.frametets.showdomain -from 0 -to 50 -increment 1 -width 3 \
            -textvariable viewoptions.drawdomainsurf -validate focus \
            -command { Ng_SetVisParameters; Ng_Vis_Set parameters; redraw } \
            -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"
            
	#tixControl $f.showdomain -label "Show surface\rof domain" -integer true \
            -variable viewoptions.drawdomainsurf -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; Ng_Vis_Set parameters; redraw }


        #tixControl $f.subdiv -label "Subdivision" -integer true \        
        #    -variable visoptions.subdivisions -min 0 -max 8 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; Ng_Vis_Set parameters; Ng_SetNextTimeStamp; redraw }
            
	#tixControl $f.frametets.showtetsdomain -label "" -integer true \
	    -variable viewoptions.drawtetsdomain -min 0 -max 500 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; redraw }

	#pack $f.frametets	
        grid $f.frametets.showtets $f.frametets.label $f.frametets.showtetsdomain -sticky w
        grid x $f.frametets.label2 $f.frametets.showdomain -stick w 
        grid x $f.frametets.label1 $f.frametets.subdiv -sticky w
        grid $f.showfilledtrigs $f.showoutline -sticky nw
        grid $f.showedges $f.showbadels -sticky nw
        grid $f.showpointnumbers $f.showedgenumbers -sticky nw
        grid $f.showfacenumbers $f.showelementnumbers -sticky nw        
        grid $f.showmetispartition $f.showidentified -sticky nw
	grid $f.showcolor $f.showpyramids -sticky nw 
        grid $f.showprisms $f.showhexes -sticky nw        
        
        
        grid  $f.frametets -sticky n -columnspan 2 -column 0 -pady 8
        #grid  $f.showdomain -stick ne;# -columnspan 3 -column 0 -pady 6
        #grid  $f.framesubdiv -sticky nw;# -columnspan 3 -column 0 -pady 6
        grid anchor $f center

        set f $w.nb.mesh
	ttk::labelframe $f.fshrink -text "Element visualization" -relief groove -borderwidth 3
	ttk::label $f.fshrink.lab -text "Shrink elements"
	#scale $f.fshrink.scale -orient horizontal -length 200 -from 0 -to 1.0001 \
	    -resolution 0.01  -tickinterval 0.25 \
	    -command { Ng_SetVisParameters; after idle redraw } \
            -variable  viewoptions.shrink
	ttk::scale $f.fshrink.scale -orient horizontal -length 200 -from 0 -to 1.0001 \
            -command "roundscale $f.fshrink.scale 2;Ng_SetVisParameters; after idle redraw" \
            -variable  viewoptions.shrink
        ttk::entry $f.fshrink.entry -textvariable viewoptions.shrink -width 4 -validate focus \
            -takefocus 0 -validatecommand "Ng_SetVisParameters; after idle redraw;my_validate %W [$f.fshrink.scale cget -from] [$f.fshrink.scale cget -to] %P 2" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; after idle redraw;"
	pack $f.fshrink -fill x -ipady 8
	grid $f.fshrink.scale $f.fshrink.entry $f.fshrink.lab -padx 4
        grid anchor $f.fshrink center
	
#	if {$userlevel == 3} {
#	    frame $f.framecurveproj
#	    checkbutton $f.framecurveproj.showcurveproj -text "Show curved edge projection " \
		-variable viewoptions.drawcurveproj \
		-command { Ng_SetVisParameters; redraw }
	#    tixControl $f.framecurveproj.showcurveprojedge -label "" -integer true \
		-variable viewoptions.drawcurveprojedge -min 1 -max 99999 \
		-options { entry.width 5 } \
		-command { Ng_SetVisParameters; redraw }

#	    pack $f.framecurveproj
#	    pack $f.framecurveproj.showcurveproj $f.framecurveproj.showcurveprojedge -side left
#	}
	
	




	# light options
	set f $w.nb.light
	ttk::labelframe $f.main -text "Lighting options" -relief groove -borderwidth 3
        pack $f.main -fill x -pady 15
        set f $f.main
	ttk::label $f.lab1 -text "Ambient Light"
	ttk::scale $f.scale1 -orient horizontal -length 200 -from 0 -to 1 \
	    -command "roundscale $f.scale1 2; Ng_SetVisParameters; redraw" \
            -variable viewoptions.light.amb
        ttk::entry $f.ent1 -textvariable viewoptions.light.amb -validate focus -width 4 \
            -validatecommand " Ng_SetVisParameters; redraw;my_validate %W [$f.scale1 cget -from] [$f.scale1 cget -to] %P 2" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"

        ttk::label $f.lab2 -text "Diffuse Light"
        ttk::scale $f.scale2 -orient horizontal -length 200 -from 0 -to 1 \
	    -command "roundscale $f.scale2 2; Ng_SetVisParameters; redraw " \
            -variable  viewoptions.light.diff 
        ttk::entry $f.ent2 -textvariable viewoptions.light.diff -validate focus -width 4 \
            -validatecommand " Ng_SetVisParameters; redraw;my_validate %W [$f.scale2 cget -from] [$f.scale2 cget -to] %P 2" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"
            
	ttk::label $f.lab3 -text "Specular Light"
	ttk::scale $f.scale3 -orient horizontal -length 200 -from 0 -to 1 \
	    -command "roundscale $f.scale3 2; Ng_SetVisParameters; redraw " \
            -variable  viewoptions.light.spec 
        ttk::entry $f.ent3 -textvariable viewoptions.light.spec -validate focus -width 4 \
            -validatecommand " Ng_SetVisParameters; redraw;my_validate %W [$f.scale3 cget -from] [$f.scale3 cget -to] %P 2" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"
        
        grid $f.scale1 $f.ent1 $f.lab1 -sticky nw -padx 4 -pady 8
        grid $f.scale2 $f.ent2 $f.lab2 -sticky nw -padx 4 -pady 8    
        grid $f.scale3 $f.ent3 $f.lab3 -sticky nw -padx 4 -pady 8
        grid anchor $f center
        set f $w.nb.light
        ttk::labelframe $f.main1 -text "Material options" -relief groove -borderwidth 3
        pack $f.main1 -fill x -pady 15
        set f $f.main1
        ttk::label $f.lab4 -text "Material Shininess"        
	ttk::scale $f.scale4 -orient horizontal -length 200 -from 0 -to 128 \
	    -command "roundscale $f.scale4 0; Ng_SetVisParameters; redraw " \
            -variable  viewoptions.mat.shininess 
        ttk::entry $f.ent4 -textvariable viewoptions.mat.shininess -validate focus -width 4 \
            -validatecommand " Ng_SetVisParameters; redraw;my_validate %W [$f.scale4 cget -from] [$f.scale4 cget -to] %P 0" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"

            
	ttk::label $f.lab5 -text "Material Transparency"
	ttk::scale $f.scale5 -orient horizontal -length 200 -from 0 -to 1 \
	-command "roundscale $f.scale5 2; Ng_SetVisParameters; redraw " \
        -variable  viewoptions.mat.transp 
        ttk::entry $f.ent5 -textvariable viewoptions.mat.transp -validate focus -width 4 \
            -validatecommand " Ng_SetVisParameters; redraw;my_validate %W [$f.scale5 cget -from] [$f.scale5 cget -to] %P 2" \
            -invalidcommand "my_invalid %W;Ng_SetVisParameters; redraw;"


        grid $f.scale4 $f.ent4 $f.lab4 -sticky nw -padx 4 -pady 8          
        grid $f.scale5 $f.ent5 $f.lab5 -sticky nw -padx 4 -pady 8
        grid anchor $f center
        #$f.lab2 $f.scale2 $f.lab3 $f.scale3 $f.lab4 $f.scale4 $f.lab5 $f.scale5
	




	# edges options
	set f $w.nb.edges
        ttk::labelframe $f.main -text "Edge viewing options" -relief groove -borderwidth 3
        pack $f.main -fill x -pady 15
        set f $f.main
        ttk::frame $f.helper 
        pack $f.helper -anchor center
        set f $f.helper
	ttk::checkbutton $f.showedges -text "Show Edges" \
	    -variable viewoptions.drawededges \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showpoints -text "Show Points" \
	    -variable viewoptions.drawedpoints \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showpointnrs -text "Show Points Nrs" \
	    -variable viewoptions.drawedpointnrs \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.showtang -text "Show CP Tangents" \
	    -variable viewoptions.drawedtangents \
	    -command { Ng_SetVisParameters; redraw }
	ttk::checkbutton $f.drawedgenrs -text "Show Edge Nrs" \
	    -variable viewoptions.drawededgenrs \
	    -command { Ng_SetVisParameters; redraw }
	
	pack $f.showedges $f.showpoints $f.showpointnrs $f.showtang $f.drawedgenrs -anchor w
        set f $w.nb.edges
        ttk::labelframe $f.main1 -text "Center point" -relief groove -borderwidth 3
        pack $f.main1 -fill x -pady 15
        set f $f.main1
	ttk::frame $f.center
	pack $f.center -anchor center 
	ttk::button $f.center.btn -text "Set Center Point" \
	    -command { Ng_SetVisParameters; Ng_Center; redraw }
	ttk::entry $f.center.ent -width 5 -textvariable viewoptions.centerpoint -validate focus \
            -validatecommand "my_validate %W 0 1e9 %P 0" \
            -invalidcommand "my_invalid %W"
	grid $f.center.ent $f.center.btn -sticky nw -padx 4
	


	#ttk::frame $f.f1
	#pack $f.f1 -pady 5 -anchor center
	ttk::label $f.center.lab1 -text "SpecPoint Veclen"
	ttk::entry $f.center.ent1 -width 5 -textvariable viewoptions.specpointvlen -validate focus \
            -validatecommand "my_validate %W 0 1e9 %P 1" \
            -invalidcommand "my_invalid %W"
	grid $f.center.ent1 $f.center.lab1 -sticky nw -padx 4
	



	# misc options
	set f $w.nb.misc

	ttk::labelframe $f.point -relief groove -borderwidth 3 -text "Special point"

	ttk::frame $f.point.dp
	
	ttk::checkbutton $f.point.dp.drawpoint -text "Draw Point" \
	    -variable viewoptions.drawspecpoint \
	    -command { Ng_SetVisParameters; redraw }

	ttk::entry $f.point.dp.px -width 8 -textvariable viewoptions.specpointx -validate focus \
            -validatecommand "my_validate %W -1e9 1e9 %P 10" \
            -invalidcommand "my_invalid %W"
	ttk::entry $f.point.dp.py -width 8 -textvariable viewoptions.specpointy -validate focus \
            -validatecommand "my_validate %W -1e9 1e9 %P 10" \
            -invalidcommand "my_invalid %W"
	ttk::entry $f.point.dp.pz -width 8 -textvariable viewoptions.specpointz -validate focus \
            -validatecommand "my_validate %W -1e9 1e9 %P 10" \
            -invalidcommand "my_invalid %W"

	grid $f.point.dp.drawpoint $f.point.dp.px $f.point.dp.py $f.point.dp.pz -sticky nw -padx 4;# -side left

	

	ttk::checkbutton $f.point.dp.center -text "Use as Center" \
	    -variable viewoptions.usecentercoords \
	    -command { 
		if { ${viewoptions.usecentercoords} } {
		    set viewoptions.centerx ${viewoptions.specpointx}
		    set viewoptions.centery ${viewoptions.specpointy}
		    set viewoptions.centerz ${viewoptions.specpointz}
		    Ng_SetVisParameters; Ng_Center
		    redraw
		} {
		    Ng_SetVisParameters
		}
		
		    
	    }

	grid $f.point.dp.center -sticky nw -padx 4
	pack $f.point.dp
	pack $f.point -fill x -ipady 3 -pady 15


	
	ttk::frame $w.bu
	pack $w.bu -fill x -ipady 3


	ttk::button $w.bu.done -text "Done" -command {
	    Ng_SetVisParameters;
	    redraw
	    destroy .viewopts_dlg
	}
	ttk::button $w.bu.apply -text "Apply" -command {
	    Ng_SetVisParameters;
	    redraw
	}
	pack $w.bu.apply $w.bu.done -expand yes -side left
	
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Viewing options"
	focus $w
    }
}



proc clipplanecommand { { optionalvar 0 } } {
    Ng_SetVisParameters
    after idle redraw 
}


set clippingdialog_pop1 0
set clippingdialog_pop2 0
set clippingdialog_pop3 0
set clippingdialog_pop4 0


#
#
#  clipping dialog
#
# 
proc clippingdialog { } {

    global clippingdialog_pop1
    global clippingdialog_pop2
    global clippingdialog_pop3
    global clippingdialog_pop4
    set clippingdialog_pop1 1
    set clippingdialog_pop2 1
    set clippingdialog_pop3 1
    set clippingdialog_pop4 1
    
    set w .clipping_dlg
    
    if {[winfo exists .clipping_dlg] == 1} {

	wm withdraw $w
	wm deiconify $w
	focus $w 

    } {
	toplevel $w
        ttk::frame $w.background
        pack $w.background -fill x -fill y 
        set w $w.background
        ttk::labelframe $w.main -text "Visual clipping" -relief groove -borderwidth 3
        pack $w.main -fill x -pady 15
        set w $w.main
	ttk::label $w.lab1 -text "Normal x"
	ttk::scale $w.scale1 -orient horizontal -length 300 -from -1 -to 1 \
            -variable  viewoptions.clipping.nx \
            -command "roundscale $w.scale1 2; clipplanecommand "
        ttk::entry $w.entry1 -width 5 -textvariable  viewoptions.clipping.nx \
            -validate focus -validatecommand " clipplanecommand;my_validate %W [$w.scale1 cget -from] [$w.scale1 cget -to] %P 2" \
            -invalidcommand "my_invalid %W; clipplanecommand"
	
	ttk::label $w.lab2 -text "Normal y"
	ttk::scale $w.scale2 -orient horizontal -length 300 -from -1 -to 1 \
            -variable  viewoptions.clipping.ny \
            -command "roundscale $w.scale2 2; clipplanecommand "
        ttk::entry $w.entry2 -width 5 -textvariable  viewoptions.clipping.ny \
            -validate focus -validatecommand " clipplanecommand;my_validate %W [$w.scale2 cget -from] [$w.scale2 cget -to] %P 2" \
            -invalidcommand "my_invalid $w.entry2;clipplanecommand"            

	ttk::label $w.lab3 -text "Normal z"
	ttk::scale $w.scale3 -orient horizontal -length 300 -from -1 -to 1 \
            -variable  viewoptions.clipping.nz \
            -command "roundscale $w.scale3 2; clipplanecommand "
        ttk::entry $w.entry3 -width 5 -textvariable  viewoptions.clipping.nz \
            -validate focus -validatecommand " clipplanecommand;my_validate %W [$w.scale3 cget -from] [$w.scale3 cget -to] %P 2" \
            -invalidcommand "my_invalid %W;clipplanecommand"

        ttk::label $w.lab4 -text "Distance"
	ttk::scale $w.scale4 -orient horizontal -length 300 -from -1 -to 1.001 \
            -variable  viewoptions.clipping.dist \
            -command "roundscale $w.scale4 3; clipplanecommand "
        ttk::entry $w.entry4 -width 5 -textvariable  viewoptions.clipping.dist \
            -validate focus -validatecommand " clipplanecommand;my_validate %W [$w.scale4 cget -from] [$w.scale4 cget -to] %P 3" \
            -invalidcommand "my_invalid %W;clipplanecommand"

        
        proc my_Press {w x y} {
            set inc [expr {([$w get $x $y] <= [$w get]) ? -1 : 1}]
            ttk::Repeatedly ttk::scale::Increment $w [expr 0.001*$inc]
            
        }
        bind $w.scale4 <ButtonPress-1> { if { [string match *slider [%W identify %x %y]] == 0 } { my_Press %W %x %y;break } }
        bind $w.scale4 <ButtonRelease-1> {ttk::scale::Release %W %x %y}
            
	ttk::label $w.lab5 -text "Additional\rDistance"
	ttk::scale $w.scale5 -orient horizontal -length 300 -from -1 -to 1.001 \
            -variable  viewoptions.clipping.dist2 \
            -command "roundscale $w.scale5 3; clipplanecommand "
        ttk::entry $w.entry5 -width 5 -textvariable  viewoptions.clipping.dist2 \
            -validate focus -validatecommand " clipplanecommand;my_validate %W [$w.scale5 cget -from] [$w.scale5 cget -to] %P 3" \
            -invalidcommand "my_invalid %W;clipplanecommand"

        bind $w.scale5 <ButtonPress-1> { if { [string match *slider [%W identify %x %y]] == 0 } { my_Press %W %x %y;break } }
        bind $w.scale5 <ButtonRelease-1> {ttk::scale::Release %W %x %y}
        
        ttk::label $w.clipdomainlabel -text "Clip only domain"
        ttk::spinbox $w.clipdomainspinb -from 0 -to 500 -increment 1 -width 3 \
            -textvariable viewoptions.clipping.onlydomain -validate focus \
            -command {clipplanecommand;} \
            -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        ttk::label $w.donotclipdomainlabel -text "Do not clip domain"
        ttk::spinbox $w.donotclipdomainspinb -from 0 -to 500 -increment 1 -width 3 \
            -textvariable viewoptions.clipping.notdomain -validate focus \
            -command "clipplanecommand" \
            -validatecommand "my_validatespinbox %W %P 0" \
            -invalidcommand "my_invalidspinbox %W"

        
	#tixControl $w.clipdomain -label "Clip only domain" -integer true \
	    -variable viewoptions.clipping.onlydomain -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { clipplanecommand; }
#	    -command { Ng_SetVisParameters; redraw }
	#tixControl $w.donotclipdomain -label "Do not clip domain" -integer true \
	    -variable viewoptions.clipping.notdomain -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { clipplanecommand; }
#	    -command { Ng_SetVisParameters; redraw }

	grid $w.scale1 $w.entry1 $w.lab1 -sticky nw -padx 4 -pady 14
        grid $w.scale2 $w.entry2 $w.lab2 -sticky nw -padx 4 -pady 14
        grid $w.scale3 $w.entry3 $w.lab3 -sticky nw -padx 4 -pady 14
        grid $w.scale4 $w.entry4 $w.lab4 -sticky nw -padx 4 -pady 14
        grid $w.scale5 $w.entry5 $w.lab5 -sticky w -padx 4 -pady 14        
        grid $w.clipdomainlabel -sticky ne -padx 4 -pady 14
        grid $w.clipdomainspinb -sticky nw -padx 4 -pady 14 -column 1 -row 5
        grid $w.donotclipdomainlabel -sticky ne -padx 4 -pady 14
        grid $w.donotclipdomainspinb -sticky nw -padx 4 -pady 14 -column 1 -row 6
        grid anchor $w center
        #pack $w.lab2 $w.scale2 $w.lab3 $w.scale3 $w.lab4 $w.scale4 $w.lab5 $w.scale5 $w.clipdomain $w.donotclipdomain

	set w .clipping_dlg.background.main
	ttk::checkbutton $w.cb1 -text "Enable clipping" \
	    -variable viewoptions.clipping.enable \
	    -command { Ng_SetVisParameters; redraw } 
            
	grid  $w.cb1 -columnspan 2 -sticky ne 
	

	
	ttk::frame $w.bu
#	pack $w.bu -fill x
	grid $w.bu;# -fill x -ipady 3
        
	ttk::button $w.cancle -text "Done" -command "destroy .clipping_dlg"
	grid $w.cancle -columnspan 3 -pady 16 
	
	set w .clipping_dlg
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Clipping Plane"
	#    grab $w
	focus $w

#	$w.scale1 configure -command { puts "call1b"; Ng_SetVisParameters; redraw } 
#	puts "after"

	clipplanecommand
    }
}





#
#  refinement dialog
#
#
proc refinementdialog { } {

    set w .refinement_dlg
    
    if {[winfo exists .refinement_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {

	toplevel $w

	#ttk::labelframe $w.main -text "Refinement options" -relief groove -borderwidth 3
        #pack $w.main -fill x -pady 15
        #set w $w.main
	# tixControl $w.meshsize -label "max mesh-size: " -integer false \
	    # -variable options.meshize -min 1e-6 -max 1e6 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }	

	# pack $w.meshsize -anchor e

	global localh
	set localh 1
	# tixControl $w.loch -label "local mesh-size: " -integer false \
	    # -variable localh -min 1e-6 -max 1e6 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }	
	
	# pack $w.loch -anchor e

	ttk::frame $w.meshsize
    
    ttk::label $w.meshsize.l1 -text "max mesh-size: "
    ttk::spinbox $w.meshsize.sp1 -from 1e-6 -to 1e6 -textvariable options.meshsize -validate focus -validatecommand "my_validatespinbox %W %P 4" \
	    -invalidcommand "my_invalidspinbox %W" -width 6 -increment 0.1
    #pack $w.meshsize.l1 $w.meshsize.sp1  -fill x -side left

	ttk::frame $w.meshsizeloc
    #pack $w.meshsize -anchor e
    #pack $w.meshsizeloc -anchor e
    ttk::label $w.meshsizeloc.l1 -text "local mesh-size: "
    ttk::spinbox $w.meshsizeloc.sp1 -from 1e-6 -to 1e6 -textvariable localh -validate focus -validatecommand "my_validatespinbox %W %P 4" \
	    -invalidcommand "my_invalidspinbox %W" -width 6 -increment 0.1
    #pack $w.meshsizeloc.l1 $w.meshsizeloc.sp1  -expand yes -fill x
    pack $w.meshsize 
    pack $w.meshsizeloc 
    grid $w.meshsize.l1 $w.meshsize.sp1
    grid $w.meshsizeloc.l1 $w.meshsizeloc.sp1
	
	
	ttk::button $w.restface -text "Restrict H at face"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH face $localh
	    }
	ttk::button $w.restedge -text "Restrict H at edge"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH edge $localh
	    }
	ttk::button $w.restelement -text "Restrict H at element"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH element $localh
	    }
	ttk::button $w.restpoint -text "Restrict H at point"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH point $localh
	    }


	pack $w.restface $w.restedge $w.restelement $w.restpoint 



	ttk::button $w.anisoedge -text "Declare Anisotropic edge"  \
	    -command {
		Ng_Anisotropy edge 
	    }
	pack $w.anisoedge
	

	frame $w.bu
	pack $w.bu -fill x -ipady 3


	ttk::button $w.bu.cancle -text "Done" -command "destroy .refinement_dlg"
	ttk::button $w.bu.refine -text "Refine"  \
	    -command { 
#		Ng_BisectCopyMesh; 
		set oldnp 0; set newnp $status_np; 
		while { $oldnp < $newnp } {
		    set level [expr $level+1]
		    Ng_Bisect; 
		    Ng_HighOrder ${options.elementorder}
		    Ng_ReadStatus;
		redraw; 
		    set oldnp $newnp
		    set newnp $status_np
		    puts "oldnp $oldnp newnp $newnp"
		}
	    }	 
	ttk::button $w.bu.zrefine -text "Z-Refine"  \
	    -command { Ng_ZRefinement; Ng_ReadStatus; redraw; }
   
	pack $w.bu.zrefine $w.bu.refine $w.bu.cancle  -expand yes -side left
		
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Select Refinement"
	focus $w
    }
}




#
#  boundcondessing dialog
#
#
proc bcpropdialog { } {

    set w .bcprop_dlg
    
    if {[winfo exists .bcprop_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
    } {
	toplevel $w
	
	frame $w.face  -borderwidth 3
	pack $w.face -fill x
	label $w.face.lab -text "face index:"
	label $w.face.ent -text 1 -padx 4
	button $w.face.next -text "next" -command {
	    set w .bcprop_dlg;	
	    set facenr [$w.face.ent cget -text]
	    if {$facenr == [Ng_BCProp getnfd]} {
		set facenr 1 
	    } {
		set facenr [expr $facenr + 1]
	    }
	    $w.face.ent configure -text $facenr
	    Ng_BCProp setactive $facenr
	    set bcnr [Ng_BCProp getbc $facenr]
	    $w.bc.ent delete 0 end
	    $w.bc.ent insert 0 $bcnr

	    redraw
	} 
	button $w.face.prev -text "prev" -command {
	    set w .bcprop_dlg;	
	    set facenr [$w.face.ent cget -text]
	    if {$facenr == 1} {
		set facenr [Ng_BCProp getnfd]
	    } {
		set facenr [expr $facenr - 1]
	    }
	    $w.face.ent configure -text $facenr
	    Ng_BCProp setactive $facenr
	    set bcnr [Ng_BCProp getbc $facenr]
	    $w.bc.ent delete 0 end
	    $w.bc.ent insert 0 $bcnr

	    redraw
	} 
	
	
	pack $w.face.lab $w.face.ent $w.face.prev $w.face.next  -side left  
	
	frame $w.bc  -borderwidth 3
	pack $w.bc -fill x
	label $w.bc.lab -text "bc property:"
	entry $w.bc.ent -width 5 -relief sunken 
	button $w.bc.but -text "change" -command { 
	    set w .bcprop_dlg;	    
	    Ng_BCProp setbc [$w.face.ent cget -text] [$w.bc.ent get]; 
	}
	button $w.bc.but2 -text "all" -command { 
	    set w .bcprop_dlg;	    
	    Ng_BCProp setall [$w.bc.ent get]; 
	}
	pack $w.bc.lab $w.bc.ent $w.bc.but $w.bc.but2 -side left  -expand yes

	frame $w.bcname  -borderwidth 3
	pack $w.bcname -fill x
	label $w.bcname.lab -text "bc name:"
	label $w.bcname.ent -text "-"
	pack $w.bcname.lab $w.bcname.ent -side left  -expand yes
	

	frame $w.bu
	pack $w.bu -fill x -ipady 3

	button $w.bu.close -text "Close" -command { destroy .bcprop_dlg }

	pack $w.bu.close  -expand yes -side left
		
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Boundary Conditions"
    }

    focus $w 

    set facenr [Ng_BCProp getactive]
    $w.face.ent configure -text $facenr
    
    set bcnr [Ng_BCProp getbc $facenr]
    $w.bc.ent delete 0 end
    $w.bc.ent insert 0 $bcnr

    set bcname [Ng_BCProp getbcname $facenr]
    $w.bcname.ent configure -text $bcname

}




#
# Philippose - 25/07/2010
# Display the face colours currently 
# available in the mesh
#
proc currmeshcoloursdialog { } {

    set w .currmeshcolours_dlg
    
    if {[winfo exists .currmeshcolours_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

	global facecolslist

	frame $w.facecols -borderwidth 3
    
	listbox $w.facecols.list -yscroll "$w.facecols.scroll set" -selectmode single -setgrid 1 -width 32 -height 12
	scrollbar $w.facecols.scroll -command "$w.facecols.list yview"
	pack $w.facecols.scroll -side right -fill y
	pack $w.facecols.list -side left -expand yes -fill both
    
	Ng_CurrentFaceColours getcolours facecolslist
	set i 1
	foreach el $facecolslist {
	    set hel [format "%d: (%.4f %.4f %.4f)" $i [ lindex $el 0 ] [ lindex $el 1 ] [ lindex $el 2 ]]
	    incr i
	    $w.facecols.list insert end $hel }

	frame $w.bu1 -borderwidth 3
    button $w.bu1.showonly -text "show only" -command {
        Ng_CurrentFaceColours showonly [.currmeshcolours_dlg.facecols.list curselection]
        redraw
    }
    button $w.bu1.hideonly -text "hide only" -command {
        Ng_CurrentFaceColours hideonly [.currmeshcolours_dlg.facecols.list curselection]
        redraw
    }
    button $w.bu1.showalso -text "show" -command {
        Ng_CurrentFaceColours showalso [.currmeshcolours_dlg.facecols.list curselection]
        redraw
    }
    button $w.bu1.hidealso -text "hide" -command {
        Ng_CurrentFaceColours hidealso [.currmeshcolours_dlg.facecols.list curselection]
        redraw
    }
    pack $w.bu1.showonly $w.bu1.hideonly $w.bu1.showalso $w.bu1.hidealso -expand yes -fill x -padx 2 -pady 2 -side left    
    
    frame $w.bu2
    button $w.bu2.showall -text "show all" -command {
        Ng_CurrentFaceColours showall
        redraw
    }
    button $w.bu2.hideall -text "hide all" -command {
        Ng_CurrentFaceColours hideall
        redraw
    }
    pack $w.bu2.showall $w.bu2.hideall -expand yes -fill x -padx 2 -pady 2 -side left 
    
    frame $w.bu3
	button $w.bu3.close -text "close" -command {
	    destroy .currmeshcolours_dlg
	}
    pack $w.bu3.close -expand yes -fill x -pady 3 -side right


	pack $w.facecols -side top -expand yes -fill x -fill y
    pack $w.bu3 -side bottom
    pack $w.bu2 -side bottom    
    pack $w.bu1 -expand yes -fill x -side left    
    
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Inspect Mesh Colours"	
	focus $w
    }
}




#
#  Philippose - 30/01/2009
#  Local Surface Mesh Size Selection
#  (Currently only supports OCC Geometry)
#
#
proc surfacemeshsizedialog { } {

    set w .surfacemeshsize_dlg

    if {[winfo exists .surfacemeshsize_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
    } {
	toplevel $w

	frame $w.face  -borderwidth 3
	pack $w.face -fill x -padx 5
	label $w.face.lab -text "face index:"
	label $w.face.ent -text 1 -padx 4
	button $w.face.next -text "next" -command {
	    set w .surfacemeshsize_dlg;	
	    set facenr [$w.face.ent cget -text]
	    if {$facenr == [Ng_SurfaceMeshSize getnfd]} {
		set facenr 1 
	    } {
		set facenr [expr $facenr + 1]
	    }
	    $w.face.ent configure -text $facenr
	    Ng_SurfaceMeshSize setactive $facenr
	    set surfms [Ng_SurfaceMeshSize getsurfms $facenr]
	    $w.sms.ent delete 0 end
	    $w.sms.ent insert 0 $surfms

	    redraw
	} 
	button $w.face.prev -text "prev" -command {
	    set w .surfacemeshsize_dlg;
	    set facenr [$w.face.ent cget -text]
	    if {$facenr == 1} {
		set facenr [Ng_SurfaceMeshSize getnfd]
	    } {
		set facenr [expr $facenr - 1]
	    }
	    $w.face.ent configure -text $facenr
	    Ng_SurfaceMeshSize setactive $facenr
	    set surfms [Ng_SurfaceMeshSize getsurfms $facenr]
	    $w.sms.ent delete 0 end
	    $w.sms.ent insert 0 $surfms

	    redraw
	} 


	pack $w.face.lab $w.face.ent $w.face.prev $w.face.next  -side left  

	frame $w.sms  -borderwidth 3
	pack $w.sms -fill x
	label $w.sms.lab -text "max mesh size:"
	entry $w.sms.ent -width 8 -relief sunken 
	button $w.sms.but -text "change" -command { 
	    set w .surfacemeshsize_dlg;	    
	    Ng_SurfaceMeshSize setsurfms [$w.face.ent cget -text] [$w.sms.ent get];
	}
	button $w.sms.but2 -text "all" -command {
	    set w .surfacemeshsize_dlg;	    
	    Ng_SurfaceMeshSize setall [$w.sms.ent get];
	}
	pack $w.sms.lab $w.sms.ent $w.sms.but $w.sms.but2 -side left -padx 5 -expand yes

	frame $w.bu
	pack $w.bu -fill x -ipady 3

	button $w.bu.close -text "Close" -command { destroy .surfacemeshsize_dlg }

	pack $w.bu.close  -expand yes -side left

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Edit Surface Mesh Size"
    }

    focus $w 

    set facenr [Ng_SurfaceMeshSize getactive]
    $w.face.ent configure -text $facenr

    set surfms [Ng_SurfaceMeshSize getsurfms $facenr]
    $w.sms.ent delete 0 end
    $w.sms.ent insert 0 $surfms

}




#
#  METIS dialog
#
#
proc METISdialog { } {

    set w .metis_dlg
    set w.parts 64
    
    if {[winfo exists .metis_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
    } {
	toplevel $w
	
	frame $w.a -borderwidth 0
	frame $w.b -borderwidth 0
	pack $w.a $w.b

	label $w.a.lab -text "Number of partitions:"
	entry $w.a.ent -textvariable w.parts -width 4 -relief sunken

	button $w.b.start -text "Start METIS" -command { 
	    Ng_Metis ${w.parts}
	    redraw
	}
	button $w.b.cancel -text "Cancel" -command { destroy .metis_dlg }
	pack $w.a.lab $w.a.ent -side left  -expand yes
	pack $w.b.start $w.b.cancel -side left


	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "METIS Partitioning"
	focus $w
 
    }
}



#
#  STL dialog
#
proc stloptionsdialog { } {

    set w .stlopts_dlg
    
    if {[winfo exists .stlopts_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

        pack [ttk::notebook $w.nb]  -fill both -fill both -side top
	# tixNoteBook $w.nb -ipadx 6 -ipady 6
	#	$w config -bg gray
	#	$w.nb subwidget nbframe config -backpagecolor gray
	
	# Create the two tabs on the notebook. The -underline option
	# puts a underline on the first character of the labels of the tabs.
	# Keyboard accelerators will be defined automatically according
	# to the underlined character.	
	#
	
# 	$w.nb add chartopt -label "Chart Options" -underline 0
# 	#$w.nb add meshsize   -label "Mesh Size"   -underline 0
# 	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	


# 	set f [$w.nb subwidget chartopt]


# 	label $f.lab1 -text "Yellow Edges Angle ()"
# 	scale $f.scale1 -orient horizontal -length 300 \
# 	    -from 0 -to 90 -resolution 1  -tickinterval 10 \
# 	    -variable  stloptions.yangle 

# 	pack $f.lab1 $f.scale1

# 	label $f.lab2e -text "Edge Corner Angle ()"
# 	scale $f.scale2e -orient horizontal -length 360 -from 0 -to 180 \
# 	    -resolution 1  -tickinterval 20 \
# 	    -variable  stloptions.edgecornerangle 
# 	pack $f.lab2e $f.scale2e
	
# 	label $f.lab2 -text "Chart Angle ()"
# 	scale $f.scale2 -orient horizontal -length 360 -from 0 -to 180 \
# 	    -resolution 1  -tickinterval 20 \
# 	    -variable  stloptions.chartangle 
# 	pack $f.lab2 $f.scale2
	
# 	label $f.lab2b -text "Outer Chart Angle ()"
# 	scale $f.scale2b -orient horizontal -length 360 -from 0 -to 180 \
# 	    -resolution 1  -tickinterval 20 \
# 	    -variable  stloptions.outerchartangle 
# 	pack $f.lab2b $f.scale2b
	
#	frame $f.r4
#	pack $f.r4 -anchor w
#	scale $f.r4.sc -orient horizontal -length 200 -from 0.1 -to 10 \
#	    -resolution 0.1 -variable stloptions.resthatlasfac
#	checkbutton $f.r4.bu -text "Restrict h for Calc Atlas (Faster)" \
#	    -variable stloptions.resthatlasenable
#	pack $f.r4.sc $f.r4.bu -side left
	

	#set f [$w.nb subwidget meshsize]

	

#    checkbutton $w.seat -text "Use Searchtrees" \
#	-variable stloptions.usesearchtree
#   pack $w.seat




	frame $w.bu
#    pack $w.bu
	pack $w.bu -fill x -ipady 3

# -fill x

    button $w.bu.apply -text "Apply" -command { redraw; Ng_GenerateMesh 1 2}
	button $w.bu.cancle -text "Done" -command { destroy .stlopts_dlg }
    pack $w.bu.cancle  $w.bu.apply  -side left -expand yes
    

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w
	wm title $w "STL Options"
#    grab $w
	focus $w
    }
}

proc stldoctordialog { } {
    Ng_STLDoctor 0 0
    set wd .stldoctor_dlg

    if {[winfo exists .stldoctor_dlg] == 1} {
    
	wm withdraw $wd
	wm deiconify $wd
	focus $wd 
    } {
	
    toplevel $wd
    pack [ttk::notebook $wd.nb] -fill both -fill both -side top
    $wd.nb add [ttk::frame $wd.nb.general] -text "General" -underline 0
    $wd.nb add [ttk::frame $wd.nb.topology] -text "Edit Topology" -underline 5
    $wd.nb add [ttk::frame $wd.nb.edges] -text "Edit Edges" -underline 5
    $wd.nb add [ttk::frame $wd.nb.normals] -text "Edit Normals" -underline 5
    $wd.nb add [ttk::frame $wd.nb.advanced] -text "Advanced" -underline 0
    
    # tixNoteBook $wd.nb -ipadx 6 -ipady 6
    # $wd.nb add general -label "General" -underline 0
    # $wd.nb add topology -label "Edit Topology"  -underline 5
    # $wd.nb add edges -label "Edit Edges"   -underline 5
    # $wd.nb add normals -label "Edit Normals"   -underline 5
    # $wd.nb add advanced -label "Advanced"   -underline 0
    # pack $wd.nb -expand yes -fill both -padx 5 -pady 5 -side top	


    # GENERAL *****************************

    set f $wd.nb.general
    ttk::frame $f.selectframe -borderwidth 0
    #ttk::frame $f.show
    #pack $f.show -fill x
    ttk::checkbutton $f.selectframe.showtrias -text "Show STL-Triangles" \
	-variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }
    #pack $f.selectframe.showtrias -anchor w
    
    ttk::checkbutton $f.selectframe.showfilledtrias -text "Show Filled Triangles" \
	-variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }
    #pack $f.show.showfilledtrias -anchor w

    set selmodevals { 0 1 2 3 4 }
    set selmodelabs(0) "triangle" 
    set selmodelabs(1) "edge" 
    set selmodelabs(2) "point" 
    set selmodelabs(3) "line" 
    set selmodelabs(4) "line cluster" 

    # tixOptionMenu $f.selmode -label "Double Click selects :" \
	# -options {
	    # label.width  19
	    # label.anchor e
	    # menubutton.width 15
	# } 

    # foreach selmodev $selmodevals {
	# $f.selmode add command $selmodev -label $selmodelabs($selmodev)
    # }
    # $f.selmode config -variable stldoctor.selectmode
    # $f.selmode config -command { Ng_STLDoctor }
    global stldoctor.selectmode
    # pack $f.selmode
    
    ttk::label  $f.selectframe.dblcsellab -text "Double Click selects : "
    ttk::menubutton $f.selectframe.dblcselbut -menu $f.selectframe.dblcselmen -text "triangle" -width 16
    menu $f.selectframe.dblcselmen  -tearoff 0
	foreach selmode { 0 1 2 3 4 } {
	    $f.selectframe.dblcselmen add command -label $selmodelabs($selmode) \
                -command "set stldoctor.selectmode $selmode ; Ng_STLDoctor ; $f.selectframe.dblcselbut configure -text \"$selmodelabs($selmode)\""
	}
    $f.selectframe.dblcselmen invoke $selmodelabs(${stldoctor.selectmode})
    pack $f.selectframe
    grid $f.selectframe.showtrias -sticky nw
    grid $f.selectframe.showfilledtrias -sticky nw
    grid $f.selectframe.dblcsellab $f.selectframe.dblcselbut -sticky nw
    
    
    
    
    ttk::frame $f.sm
    pack $f.sm -fill x
    ttk::checkbutton $f.sm.bu -text "select with mouse" \
	-variable stldoctor.selectwithmouse
    pack $f.sm.bu 

    ttk::frame $f.st -relief groove -borderwidth 3
    pack $f.st -fill x
    ttk::label $f.st.lab -text "Select triangle by number";
    ttk::entry $f.st.ent -width 5 \
	-textvariable stldoctor.selecttrig
    pack $f.st.ent $f.st.lab -side left -expand yes

    ttk::frame $f.vc -relief groove -borderwidth 3
    pack $f.vc -fill x
    ttk::checkbutton $f.vc.bu -text "show vicinity" \
	-variable stldoctor.showvicinity \
	-command {Ng_STLDoctor vicinity; redraw}
    ttk::label $f.vc.lab -text "vicinity size";
    
    #scale $f.vc.sc -orient horizontal -length 200 -from 0 -to 200 \
	-resolution 1 -variable stldoctor.vicinity \
	-command { Ng_STLDoctor vicinity; redraw }
    ttk::frame $f.vc.sc
    ttk::scale $f.vc.sc.scale -orient horizontal -length 200 -from 0 -to 200 \
    -variable stldoctor.vicinity -takefocus 0 -command "Ng_STLDoctor vicinity; redraw; roundscale $f.vc.sc.scale 0"
    ttk::entry $f.vc.sc.entry -textvariable stldoctor.vicinity -width 3 \
    -validatecommand "Ng_STLDoctor vicinity; redraw; my_validate %W [$f.vc.sc.scale cget -from] [$f.vc.sc.scale cget -to] %P 0" \
    -invalidcommand "my_invalid %W;Ng_STLDoctor vicinity; redraw;" -validate focus
    ttk::label $f.vc.sc.lab -text "vicinity size"
	grid $f.vc.sc.scale $f.vc.sc.entry $f.vc.sc.lab -sticky nw -padx 4    
    pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes

    ttk::frame $f.ge -relief groove -borderwidth 0
    pack $f.ge -expand yes
    ttk::button $f.ge.neighbourangles -text "calc neighbourangles" -command {Ng_STLDoctor neighbourangles}
    ttk::button $f.ge.showcoords -text "show coords of touched triangle" -command {Ng_STLDoctor showcoords}
    ttk::button $f.ge.moveptm -text "move point to middle of trianglepoints" -command {Ng_STLDoctor movepointtomiddle; redraw}
    ttk::button $f.ge.destroy0trigs -text "destroy 0-volume triangles" -command {Ng_STLDoctor destroy0trigs}
    grid $f.ge.neighbourangles -sticky nw -padx 4 -pady 4
    grid $f.ge.showcoords -sticky nw -padx 4 -pady 4
    grid $f.ge.moveptm -sticky nw -padx 4 -pady 4
    grid $f.ge.destroy0trigs -sticky nw -padx 4 -pady 4


    ttk::button $f.ge.cancle -text "Done" -command {destroy .stldoctor_dlg }
    grid $f.ge.cancle  -sticky nw

    # TOPOLOGY ********************
    set f $wd.nb.topology

    ttk::frame $f.oc -relief groove -borderwidth 3
    pack $f.oc -pady 3 -ipady 3 -fill y -fill x
    ttk::frame $f.oc.oc1 -borderwidth 0
    pack $f.oc.oc1
    ttk::button $f.oc.oc1.bu -text "invert orientation \n of selected trig" -command {Ng_STLDoctor invertselectedtrig; redraw }
    ttk::button $f.oc.oc1.bu2 -text "orient after \n selected trig" -command {Ng_STLDoctor orientafterselectedtrig; redraw }
    

    ttk::button $f.oc.oc1.toperr -text "mark inconsistent triangles" -command {Ng_STLDoctor marktoperrortrigs; redraw }
    ttk::button $f.oc.oc1.deltrig -text "delete selected triangle" -command {Ng_STLDoctor deleteselectedtrig; redraw }
    ttk::button $f.oc.oc1.geosmooth -text "geometric smoothing" -command {Ng_STLDoctor smoothgeometry; redraw }
    
    grid $f.oc.oc1.bu x $f.oc.oc1.bu2 -sticky nw -padx 4 -pady 4
    grid $f.oc.oc1.toperr - x -sticky nw -padx 4 -pady 4
    grid $f.oc.oc1.deltrig - x -sticky nw -padx 4 -pady 4
    grid $f.oc.oc1.geosmooth - x -sticky nw -padx 4 -pady 4





    # EDGES ***********************
    set f $wd.nb.edges


    ttk::frame $f.be -relief groove -borderwidth 3 
    pack $f.be -fill x
    
    #scale $f.be.sc -orient horizontal -length 200 -from 0 -to 100 \
	#-resolution 0.5
    ttk::frame $f.be.frame
    pack $f.be.frame -ipady 4 -pady 4
    ttk::label $f.be.frame.lab -text "build edges with yellow angle:";
    ttk::scale $f.be.frame.scale -orient horizontal -length 200 -from 0 -to 200 \
    -variable stloptions.yangle -takefocus 0 -command "roundscale $f.be.frame.scale 1; Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw"
    ttk::entry $f.be.frame.entry -textvariable stloptions.yangle -width 5 \
    -validatecommand "Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw;my_validate %W [$f.be.frame.scale cget -from] [$f.be.frame.scale cget -to] %P 1" \
    -invalidcommand "my_invalid %W;Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw" -validate focus
	grid $f.be.frame.lab - -sticky nw -padx 4
    grid $f.be.frame.scale $f.be.frame.entry -sticky nw -padx 4    
    
    
    
    #$f.be.sc config -variable stloptions.yangle 
    #$f.be.sc config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }
    ttk::label $f.be.frame.lab2 -text "continue edges with yellow angle:";
#    scale $f.be.sc2 -orient horizontal -length 200 -from 0 -to 100 \
	-resolution 0.5
    ttk::scale $f.be.frame.scale2 -orient horizontal -length 200 -from 0 -to 100 \
    -variable stloptions.contyangle -takefocus 0 -command "roundscale $f.be.frame.scale2 1; Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw"
    ttk::entry $f.be.frame.entry2 -textvariable stloptions.contyangle -width 5 \
    -validatecommand "Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw;my_validate %W [$f.be.frame.scale2 cget -from] [$f.be.frame.scale2 cget -to] %P 1" \
    -invalidcommand "my_invalid %W;Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw" -validate focus

	grid $f.be.frame.lab2 - -sticky nw -padx 4
    grid $f.be.frame.scale2 $f.be.frame.entry2 -sticky nw -padx 4    
    
    #$f.be.sc2 config -variable stloptions.contyangle 
    #$f.be.sc2 config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }



    ttk::button $f.be.frame.buildedges -text "Build Edges" -command {Ng_STLDoctor buildedges; redraw}
    grid $f.be.frame.buildedges - -sticky n -padx 4 -pady 4
    #pack $f.be.lab $f.be.sc $f.be.lab2 $f.be.sc2 $f.be.buildedges -expand yes

    ttk::frame $f.se -relief groove -borderwidth 3
    pack $f.se -fill x
    ttk::checkbutton $f.se.bu -text "show excluded" \
	-variable stldoctor.showexcluded \
	-command {Ng_STLDoctor; redraw}
    pack $f.se.bu 

    # edgeselectmode ******

    set edgeselmodevals { 0 1 2 3 4 }
    set edgeselmodelabs(0) "no change" 
    set edgeselmodelabs(1) "undefined" 
    set edgeselmodelabs(2) "confirmed" 
    set edgeselmodelabs(3) "candidate"
    set edgeselmodelabs(4) "excluded"

    # tixOptionMenu $f.edgeselmode -label "Double Click sets edge :" \
	# -options {
	    # label.width  19
	    # label.anchor e
	    # menubutton.width 15
	# } 

    # foreach edgeselmodev $edgeselmodevals {
	# $f.edgeselmode add command $edgeselmodev -label $edgeselmodelabs($edgeselmodev)
    # }
    # $f.edgeselmode config -variable stldoctor.edgeselectmode
    # $f.edgeselmode config -command { Ng_STLDoctor }
    global stldoctor.edgeselectmode    
    # pack $f.edgeselmode

    ttk::frame $f.scaleframe -relief groove -borderwidth 0
    pack $f.scaleframe -ipadx 4 -pady 4 -expand yes
    ttk::label  $f.scaleframe.dblcedgelab -text "Double Click sets edge :"
    ttk::menubutton $f.scaleframe.dblcledgebut -menu $f.scaleframe.dblcledgem -text "coarse" -width 16

    menu $f.scaleframe.dblcledgem  -tearoff 0
	foreach selectmode { 0 1 2 3 4 } {
	$f.scaleframe.dblcledgem add command -label $edgeselmodelabs($selectmode) \
    -command "set stldoctor.edgeselectmode $selectmode ; $f.scaleframe.dblcledgebut configure -text \"$edgeselmodelabs($selectmode)\""
	}
    $f.scaleframe.dblcledgem invoke $edgeselmodelabs(${stldoctor.edgeselectmode})                
    grid $f.scaleframe.dblcedgelab $f.scaleframe.dblcledgebut -sticky n -ipadx 4

    
    
    
    # edge buttons

    ttk::frame $f.edg -relief groove -borderwidth 3
    pack $f.edg -fill x -ipadx 4 -ipady 4

#    checkbutton $f.edg.bu -text "use external edges" \
#	-variable stldoctor.useexternaledges \
#	-command {Ng_STLDoctor; redraw}
#   pack $f.edg.bu -expand yes


    ttk::frame $f.edg.f0
    pack $f.edg.f0
    ttk::button $f.edg.f0.confirmedge -text "confirm" -command {Ng_STLDoctor confirmedge; redraw}
    ttk::button $f.edg.f0.candidateedge -text "candidate" -command {Ng_STLDoctor candidateedge; redraw}
    ttk::button $f.edg.f0.excludeedge -text "exclude" -command {Ng_STLDoctor excludeedge; redraw}
    ttk::button $f.edg.f0.undefinededge -text "undefined" -command {Ng_STLDoctor undefinededge; redraw}
    pack $f.edg.f0.confirmedge $f.edg.f0.candidateedge $f.edg.f0.excludeedge $f.edg.f0.undefinededge  -side left

    ttk::frame $f.edg.fa
    pack $f.edg.fa
    ttk::button $f.edg.fa.setallundefined -text "all undefined" -command {Ng_STLDoctor setallundefinededges; redraw}
    ttk::button $f.edg.fa.erasecandidates -text "candidates to undefined" -command {Ng_STLDoctor erasecandidateedges; redraw}
    pack $f.edg.fa.setallundefined $f.edg.fa.erasecandidates -side left


    ttk::frame $f.edg.fb
    pack $f.edg.fb
    ttk::button $f.edg.fb.confirmcandidates -text "candidates to confirmed" -command {Ng_STLDoctor confirmcandidateedges; redraw}
    ttk::button $f.edg.fb.confirmedtocandidates -text "confirmed to candidates" -command {Ng_STLDoctor confirmedtocandidateedges; redraw}
    pack $f.edg.fb.confirmcandidates $f.edg.fb.confirmedtocandidates -side left

    ttk::frame $f.edg.f1
    ttk::frame $f.edg.f2
    ttk::frame $f.edg.f3
    ttk::frame $f.edg.f4
    pack $f.edg.f1 $f.edg.f2 $f.edg.f3 $f.edg.f4

    ttk::button $f.edg.f1.exportedges -text "export edges" -command {Ng_STLDoctor exportedges}
    ttk::button $f.edg.f1.importedges -text "import edges" -command {Ng_STLDoctor importedges; redraw}
    ttk::button $f.edg.f1.saveedgedata -text "save edgedata" \
	-command { 
	    set types {
		{"Netgen Edgedata"   {.ned} } 
	    }
	    set file [tk_getSaveFile -filetypes $types -defaultextension ".ned"]
	    if {$file != ""} {
		Ng_STLDoctor saveedgedata $file
	}
    }

    ttk::button $f.edg.f1.loadedgedata -text "load edgedata" \
	-command { 
	    set types {
		{"Netgen Edgedata"  {.ned} }
	    }
	    set file [tk_getOpenFile -filetypes $types -defaultextension ".ned"]
	    if {$file != ""} {
		Ng_STLDoctor loadedgedata $file 
		puts "loading done"
		
		redraw
		
#		wm title . [concat "NETGEN - " $file]
	    }
	} 

    ttk::button $f.edg.f1.importAVLedges -text "import AVL edges" \
	-command {
	    set types {{"Edge file"  {.edg }}}

	    set file [tk_getOpenFile -filetypes $types -defaultextension ".edg"]
	    if {$file != ""} {
		Ng_STLDoctor importexternaledges $file; 
	    }
	}

    pack $f.edg.f1.importAVLedges $f.edg.f1.loadedgedata $f.edg.f1.saveedgedata -side left

#    button $f.edg.f1.buildedges -text "build external edges" -command {Ng_STLDoctor buildexternaledges; redraw}
    ttk::frame $f.edg2 -relief groove -borderwidth 3
    pack $f.edg2 -fill x


#    button $f.edg2.addlonglines -text "make long lines candidates (% of diam)" -command {Ng_STLDoctor addlonglines; redraw}
    label $f.edg2.lab -text "length (%):"
    scale $f.edg2.sc -orient horizontal -length 200 -from 0 -to 100 \
	-resolution 0.5 \
        -variable stldoctor.longlinefact 

 #   button $f.edg2.deletedirtyedges -text "make dirty edges candidates" -command {Ng_STLDoctor deletedirtyedges; redraw}
    button $f.edg2.undoedge -text "undo last edge change" -command {Ng_STLDoctor undoedgechange; redraw}
    
 #   pack $f.edg2.addlonglines $f.edg2.deletedirtyedges -expand yes
 #   pack $f.edg2.lab $f.edg2.sc -side left
    pack $f.edg2.undoedge -expand yes



    # NORMALS ***********************
    set f $wd.nb.normals

    ttk::frame $f.dt -relief groove -borderwidth 3
    pack $f.dt -fill x
    ttk::label $f.dt.lab -text "dirty triangle factor";
    ttk::entry $f.dt.ent -width 5 \
	-textvariable stldoctor.dirtytrigfact -validatecommand "Ng_SetSTLParameters;my_validate %W -1e9 1e9 %P 3" \
    -invalidcommand "my_invalid %W;Ng_SetSTLParameters" -validate focus
    pack $f.dt.ent $f.dt.lab -side left -expand yes -pady 8

    ttk::frame $f.srt -relief groove -borderwidth 3
    pack $f.srt -fill x
    ttk::button $f.srt.bu -text "smooth reverted triangles geometric" -command {Ng_STLDoctor smoothrevertedtrigs; redraw }
    ttk::entry $f.srt.ent -width 5 \
	-textvariable stldoctor.smoothangle -validatecommand "Ng_SetSTLParameters;my_validate %W -1e9 1e9 %P 2" \
    -invalidcommand "my_invalid %W;Ng_SetSTLParameters" -validate focus
    pack $f.srt.ent $f.srt.bu -side left  -expand yes -pady 8

    ttk::frame $f.bdt -relief groove -borderwidth 3
    pack $f.bdt -fill x
    ttk::button $f.bdt.bu -text "mark dirty triangles" -command {Ng_STLDoctor markdirtytrigs; redraw }
    ttk::button $f.bdt.bu2 -text "smooth dirty triangles normal" -command {Ng_STLDoctor smoothdirtytrigs; redraw }
    pack $f.bdt.bu $f.bdt.bu2 -side left  -expand yes -pady 8

    
    ttk::frame $f.sno -relief groove -borderwidth 3
    pack $f.sno -fill x
    ttk::frame $f.sno.snoframe -borderwidth 0
    #ttk::label $f.sno.labrough -text "rough"
    #scale $f.sno.scsmooth -orient horizontal -length 100 -from 0 -to 0.8 \
	-resolution 0.01 -variable stldoctor.smoothnormalsweight \
	-command { Ng_SetSTLParameters }
    #ttk::label $f.sno.labsmooth -text "smooth"
    ttk::button $f.sno.smoothnormals -text "smooth normals" -command { Ng_STLDoctor smoothnormals; redraw}

    ttk::scale $f.sno.snoframe.scale -orient horizontal -length 100 -from 0.0 -to 0.8 \
    -variable stldoctor.smoothnormalsweight -takefocus 0 -command "roundscale $f.sno.snoframe.scale 2;Ng_SetSTLParameters"
    ttk::entry $f.sno.snoframe.entry -textvariable stldoctor.smoothnormalsweight -width 4 \
    -validatecommand "Ng_SetSTLParameters;my_validate %W [$f.sno.snoframe.scale cget -from] [$f.sno.snoframe.scale cget -to] %P 2" \
    -invalidcommand "my_invalid %W;Ng_SetSTLParameters" -validate focus
    ttk::label $f.sno.snoframe.labrough -text "rough"
    ttk::label $f.sno.snoframe.labsmooth -text "smooth"
	grid $f.sno.snoframe.labrough $f.sno.snoframe.scale $f.sno.snoframe.labsmooth $f.sno.snoframe.entry -sticky nw -padx 4    
    
    #pack $f.sno.labrough $f.sno.scsmooth $f.sno.labsmooth $f.sno.smoothnormals -side left -padx 5
    pack $f.sno.snoframe $f.sno.smoothnormals -side left -padx 5 -pady 8
    ttk::frame $f.no -relief groove -borderwidth 3
    pack $f.no -fill x

    ttk::button $f.no.marknonsmoothnormals -text "mark non-smooth triangles" -command {Ng_STLDoctor marknonsmoothnormals; redraw}
    ttk::button $f.no.calcnormals -text "calculate normals from geometry" -command {Ng_STLDoctor calcnormals; redraw}

    pack $f.no.marknonsmoothnormals $f.no.calcnormals -expand yes -pady 8


    # ADVANCED **************************
    set f $wd.nb.advanced


    ttk::frame $f.sc
    pack $f.sc -fill x
    ttk::checkbutton $f.sc.bu -text "spiral check" \
	-variable stldoctor.spiralcheck \
	-command {Ng_STLDoctor;}    
    ttk::checkbutton $f.sc.bu2 -text "cone check" \
	-variable stldoctor.conecheck \
	-command {Ng_STLDoctor;}    
    pack $f.sc.bu $f.sc.bu2

    
    #tixControl $f.gtol -label "load-geometry tolerance factor" -integer false \
	-variable stldoctor.geom_tol_fact \
	-options {
	#    entry.width 8
	#    label.width 30
	#    label.anchor e
	#}	
    ttk::spinbox $f.gtol -from 1 -to 20 -textvariable stldoctor.geom_tol_fact -width 8
    pack $f.gtol

    ttk::button $f.adap -text "Apply" -command {
	.stldoctor_dlg.nb.advanced.gtol invoke
	Ng_STLDoctor; 
    }
    pack $f.adap -expand yes

#    frame $f.gtol -relief groove -borderwidth 3
#    pack $f.gtol -fill x
#    label $f.gtol.lab -text "Geometry-Load-Tolerance-Factor";
#    entry $f.gtol.ent -width 5 -relief sunken \
#	-textvariable stldoctor.geom_tol_fact
#   pack $f.gtol.lab $f.gtol.ent -side left -expand yes

    #*******************************    
    wm withdraw $wd
    wm geom $wd +100+100
    wm deiconify $wd
    wm title $wd "STL Doctor"

    focus $wd
}
}





proc meshdoctordialog { } {

    set w .meshdoc_dlg
    global meshdoctor.active

    if {[winfo exists .meshdoc_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

	set meshdoctor.active 1
	Ng_MeshDoctor;


	ttk::frame $w.vis -relief groove -borderwidth 3
	pack $w.vis

	ttk::checkbutton $w.vis.showfilledtrigs -text "Show filled triangles" \
	-variable viewoptions.drawfilledtrigs \
	-command { Ng_SetVisParameters; redraw }
	
	ttk::checkbutton $w.vis.showedges -text "Show edges" \
	-variable viewoptions.drawedges \
	-command { Ng_SetVisParameters; redraw }
	

	ttk::checkbutton $w.vis.showoutline -text "Show Triangle Outline" \
	-variable viewoptions.drawoutline \
	-command { Ng_SetVisParameters; redraw }

	pack $w.vis.showfilledtrigs  $w.vis.showoutline $w.vis.showedges

    ttk::frame $w.markedgedist
    ttk::label $w.markedgedist.l -text "Mark edge dist: "
    ttk::spinbox $w.markedgedist.s -from 0 -to 999 -width 5 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
    -invalidcommand "my_invalidspinbox %W" -command {Ng_MeshDoctor markedgedist ${meshdoc.markedgedist};redraw} -textvariable meshdoc.markedgedist
    #pack $f.grading -fill x
    pack $w.markedgedist.l $w.markedgedist.s -side left
    
	# tixControl $w.markedgedist -label "Mark edge dist: " -integer true \
	    # -min 0 -max 999  \
	    # -variable meshdoc.markedgedist \
	    # -options {
		# entry.width 3
		# label.width 20
		# label.anchor e
	    # } \
	    # -command {
		# Ng_MeshDoctor markedgedist ${meshdoc.markedgedist}
		# redraw
	    # }
	pack $w.markedgedist
	
	ttk::button $w.deledge -text "Delete marked segments" -command {
	    Ng_MeshDoctor deletemarkedsegments
	    redraw
	}
	pack $w.deledge
	
	ttk::button $w.close -text "Close" -command { 
	    set meshdoctor.active 0;
	    Ng_MeshDoctor;
	    destroy .meshdoc_dlg 
	}
	pack $w.close -expand yes
	
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Mesh Doctor"
    }
}



#
#  Quality viewer
#

proc qualityviewdialog { show } {

    set w .qualityview_dlg
    
    if {[winfo exists .qualityview_dlg] == 1} {

	if { $show == 1 } {
	    wm withdraw .qualityview_dlg
	    wm deiconify $w
	    focus $w 
	} {
	    wm withdraw $w
	}
    } {
	toplevel $w
	
	set c $w.c

	canvas $c -relief raised -width 450 -height 300
	pack $w.c -side top -fill x

	set plotFont {Helvetica 12}
	set smallFont {Helvetica 12}

	$c create line 100 250 400 250 -width 2
	$c create line 100 250 100 50 -width 2

	for {set i 0} {$i <= 10} {incr i} {
	    set x [expr {100 + ($i*30)}]
	    $c create line $x 250 $x 245 -width 2
	    if { [expr {$i % 2}] == 0 } {
		$c create text $x 254 -text [format %1.1f [expr 0.1*$i]] -anchor n -font $plotFont
	    }
	}

	global qualbar
	global qualbarnull
	global qualbaraxis

	for {set i 0} {$i <= 5} {incr i} {
	    set y [expr {250 - ($i*40)}]
	    $c create line 100 $y 105 $y -width 2

#            global qualbaraxis($i)
	    set qualbaraxis($i) \
		[$c create text 96 $y -text [expr $i*50].0 -anchor e -font $plotFont]
	}

	for {set i 0} {$i < 20} {incr i} {
	    set x1 [expr {100 + ($i*15) + 2}]
	    set x2 [expr {$x1+10}]
	    set y [expr {250 - 10 * $i}]
#	    global qualbar($i)
	    set qualbar($i) [$c create rectangle $x1 250 $x2 245 -fill blue]
	    set qualbarnull($i) [$c create text [expr {($x1+$x2)/2}] 245 -text 0 -anchor s -font $smallFont -fill blue]	
	}

	frame $w.bu
	pack $w.bu
	# -fill x
	
	button $w.close -text "Close" \
	    -command { 
		wm withdraw .qualityview_dlg
		set viewqualityplot 0
	    }
	pack $w.close
    
	
	if { $show == 1 } {
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "Mesh Quality"
	    focus $w
	}
    }
}










#
#  Quality viewer
#
proc memusedialog { show } {

    set w .memuse_dlg
    
    if {[winfo exists .memuse_dlg] == 1} {

	if { $show == 1 } {
	    wm withdraw .memuse_dlg
	    wm deiconify $w
	    focus $w 
	} {
	    wm withdraw $w
	}
    } {
	toplevel $w
	
	set c $w.c

	canvas $c -relief raised -width 600 -height 300
	pack $w.c -side top -fill x

	set plotFont {Helvetica 18}
	set smallFont {Helvetica 12}


	global memmark
	for {set i 0} {$i < 512} { incr i } {
	    set memmark($i) [$c create line [expr 50+$i] 50 [expr 50+$i] 70 -fill blue]
	}


	set plotFont {Helvetica 18}
	set smallFont {Helvetica 12}

	$c create text 50 90 -text "0 GB" -anchor n -font $plotFont
	$c create text 178 90 -text "1 GB" -anchor n -font $plotFont
	$c create text 306 90 -text "2 GB" -anchor n -font $plotFont
	$c create text 434 90 -text "3 GB" -anchor n -font $plotFont
	$c create text 562 90 -text "4 GB" -anchor n -font $plotFont


	frame $w.bu
	pack $w.bu
	# -fill x
	
	button $w.close -text "Close" \
	    -command { 
		wm withdraw .memuse_dlg
		set memuseplot 0
	    }
	pack $w.close
	
	if { $show == 1 } {
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "Memory Usage"
	    focus $w
	}
    }
}















#
#  STL INFO dialog
#
proc STLinfodialog { show } {

    set w .STLinfo_dlg
    
    if {[winfo exists .STLinfo_dlg] == 1} {

	if { $show == 1 } {
	    wm withdraw .STLinfo_dlg
	    wm deiconify $w
	    focus $w 
	} {
	    wm withdraw $w
	}
    } {
	toplevel $w
	
	set c $w.c

	canvas $c -relief raised -width 450 -height 300
	pack $w.c -side top -fill x

	set plotFont {Helvetica 18}
	set smallFont {Helvetica 12}

	$c create line 100 250 400 250 -width 2
	$c create line 100 250 100 50 -width 2

	frame $w.bu
	pack $w.bu
	# -fill x
	
	button $w.close -text "Close" \
	    -command { 
		wm withdraw .STLinfo_dlg
		#set STLinfoopen 0
	    }
	pack $w.close
    
	
	if { $show == 1 } {
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "STL Geometry Info"
	    focus $w
	}
    }
}








proc logwindow { } {
    set w .logwindow
    
    if {[winfo exists .logwindow] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

	text $w.edit -yscroll "$w.scrolly set" -setgrid 1 -height 12
	scrollbar $w.scrolly -command "$w.edit yview"	
	pack $w.edit -side left -fill both -expand 1
	pack $w.scrolly -side left -fill both -expand 0

	.logwindow.edit insert end "Netgen Log Window\n"

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Netgen Log"	
	focus $w
    }
}
# logwindow



# Opens a window with a table. tablevar is a list, the first entry is the title, the second the number of rows, the third the number of columns,
# then the entries follow.

proc printtable { tablevar } {
    set w newtcltable
    while {[winfo exists .$w] == 1} {set w 1$w}
    set w .$w
    toplevel $w
    for {set i 0} {$i < [lindex $tablevar 2]} { incr i } {
	frame $w.col$i
	for {set j 0} {$j < [lindex $tablevar 1]} { incr j } {
	    frame $w.col$i.row$j
	    message $w.col$i.row$j.txt -aspect 10000000 -text [lindex $tablevar [expr 3+[lindex $tablevar 2]*$j+$i]]
	    pack $w.col$i.row$j.txt
	    pack $w.col$i.row$j -side top
	}
	pack $w.col$i -side left
    }
    wm withdraw $w
    wm geom $w +200+100; wm deiconify $w
    wm title $w [lindex $tablevar 0]
    focus $w
}


set latestwarning 0


proc printwarning { textvar } {
    global latestwarning
    set latestwarning $textvar
    set w warning
    while {[winfo exists .$w] == 1} {set w 1$w}
    set w .$w
    toplevel $w
    message $w.mes -aspect 2000 -text "WARNING:\n$textvar"
    button $w.done -text "Done" -command "destroy $w"
    pack $w.mes
    pack $w.done
    wm withdraw $w
    wm deiconify $w
    wm title $w "Warning"
    focus $w
}


proc printlatestwarning { } { 
    global latestwarning 
    if {$latestwarning != 0} {printwarning $latestwarning}
}



# proc runtestdialog { } {
    # source $::ngdir/ngshell.tcl
    # set w .runtest_dlg
    
    # if {[winfo exists .runtest_dlg] == 1} {
	# wm withdraw $w
	# wm deiconify $w

	# focus $w 
    # } {
	# toplevel $w

# # in2d testing #
	# frame $w.in2dframe 
	# pack $w.in2dframe

        # set in2dlogfile ""
 	# tixLabelEntry $w.in2dframe.ent -label "in2d log-file: console if empty"  \
 	    # -labelside top \
 	    # -options {  
 		# entry.textVariable in2dlogfile
 		# entry.width 35
 		# label.width 25
 		# label.anchor w
 	    # }	
 	# button $w.in2dframe.btn -text "Browse" -command {
	    # set types { { "Log file"   {.log}	} }
	    # set in2dlogfile [tk_getOpenFile -filetypes $types -initialfile $in2dlogfile]
	# }
 	# button $w.in2dframe.test -text "Test in2d meshing" -command { ngtest in2d $in2dlogfile }

        
 	# pack $w.in2dframe.test -side left -anchor s -padx 4 -pady 4
 	# pack $w.in2dframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	# pack $w.in2dframe.btn -side left -anchor s -padx 4 -pady 4

        
# # geo testing #
	# frame $w.geoframe 
	# pack $w.geoframe

        # set geologfile "" 
 	# tixLabelEntry $w.geoframe.ent -label "geo log-file: console if empty"  \
 	    # -labelside top \
 	    # -options {  
 		# entry.textVariable geologfile
 		# entry.width 35
 		# label.width 25
 		# label.anchor w
 	    # }	
 	# button $w.geoframe.btn -text "Browse" -command {
	    # set types { { "Log file"   {.log}	} }
	    # set geologfile [tk_getOpenFile -filetypes $types -initialfile $geologfile]
	# }
 	# button $w.geoframe.test -text "Test geo meshing" -command { ngtest geo $geologfile }

        
 	# pack $w.geoframe.test -side left -anchor s -padx 4 -pady 4
 	# pack $w.geoframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	# pack $w.geoframe.btn -side left -anchor s -padx 4 -pady 4

# # stl testing #
	# frame $w.stlframe 
	# pack $w.stlframe

        # set stllogfile ""
 	# tixLabelEntry $w.stlframe.ent -label "stl log-file: console if empty"  \
 	    # -labelside top \
 	    # -options {  
 		# entry.textVariable stllogfile
 		# entry.width 35
 		# label.width 25
 		# label.anchor w
 	    # }	
 	# button $w.stlframe.btn -text "Browse" -command {
	    # set types { { "Log file"   {.log}	} }
	    # set stllogfile [tk_getOpenFile -filetypes $types -initialfile $stllogfile]
	# }
 	# button $w.stlframe.test -text "Test stl meshing" -command { ngtest stl $stllogfile }

        
 	# pack $w.stlframe.test -side left -anchor s -padx 4 -pady 4
 	# pack $w.stlframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	# pack $w.stlframe.btn -side left -anchor s -padx 4 -pady 4

# # pde testing #
	# frame $w.pdeframe 
	# pack $w.pdeframe

        # set pdelogfile ""
 	# tixLabelEntry $w.pdeframe.ent -label "pde log-file: console if empty"  \
 	    # -labelside top \
 	    # -options {  
 		# entry.textVariable pdelogfile
 		# entry.width 35
 		# label.width 25
 		# label.anchor w
 	    # }	
 	# button $w.pdeframe.btn -text "Browse" -command {
	    # set types { { "Log file"   {.log}	} }
	    # set pdelogfile [tk_getOpenFile -filetypes $types -initialfile $pdelogfile]
	# }
 	# button $w.pdeframe.test -text "Test ngsolve pde's" -command { ngtest pde $pdelogfile }

        
 	# pack $w.pdeframe.test -side left -anchor s -padx 4 -pady 4
 	# pack $w.pdeframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	# pack $w.pdeframe.btn -side left -anchor s -padx 4 -pady 4
 
	# wm title $w "Testing"
	# focus .runtest_dlg 
    # }
# }

