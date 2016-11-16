Ng_Vis_Set parameters



set viscnt 0
proc snapshottimer { } {
    after 2000 { snapshottimer }

    global viscnt
    set viscnt [expr $viscnt+1]
    set s1 0000$viscnt
#    puts $s1
    set cnt [string range $s1 [expr [string length $s1]-4] end]
    set filename "p$cnt.jpg"
#    puts "cnt = $cnt"
#    puts "filename = $filename"
#    .ndraw Ng_SnapShot pictures/$filename
}
snapshottimer



proc redrawtimer { } {
    global visoptions.autoredraw
    global visoptions.autoredrawtime

    set delay [expr int(${visoptions.autoredrawtime}*1000)]
    if { ${visoptions.autoredraw} == 1 } { redraw; }
    after $delay { redrawtimer } 
}
redrawtimer


set perstarttime [clock clicks -millisecond]
proc redrawperiodic { } {
    global visoptions.redrawperiodic
    global perstarttime
    set curtime [clock clicks -millisecond]
#    puts "redraw periodic, time = $curtime"
    Ng_Vis_Set time [expr ($curtime - $perstarttime) / 5]
    redraw
    if { ${visoptions.redrawperiodic} == 1 } { after 30 { redrawperiodic } };

}



proc addplotline { identifier datax datay plotinfo {color black}} {
    set c $identifier.c

    set xstart [lindex $plotinfo 0]
    set ystart [lindex $plotinfo 1]
    set xmin [lindex $plotinfo 2]
    set ymin [lindex $plotinfo 3]
    set unitx [lindex $plotinfo 4]
    set unity [lindex $plotinfo 5]
    
    
    
    set latestx [expr ([lindex $datax 0]-$xmin)*$unitx + $xstart]
    set latesty [expr ([lindex $datay 0]-$ymin)*$unity + $ystart]
    
    for {set i 1} {$i < [llength $datax]} {incr i} {
	set xpos [expr ([lindex $datax $i]-$xmin)*$unitx + $xstart]
	set ypos [expr ([lindex $datay $i]-$ymin)*$unity + $ystart]
	$c create line $latestx $latesty $xpos $ypos -width 1 -fill $color
	set latestx $xpos
	set latesty $ypos
    }
}



proc createlineplot { width height identifier title xmin xmax ymin ymax plotinfo} {
    set thiswidth $width
    set thisheight $height
    if { $thiswidth < 275 } { set thiswidth 275 }
    if { $thisheight < 225 } { seth thisheight 225 }

    set w $identifier

    if {[winfo exists $w] == 1} {
	
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w
	
	set c $w.c
	
	canvas $c -relief raised -width $thiswidth -height $thisheight
	pack $w.c -side top -fill x
	
	set titleFont {Helvetica 18}
	set smallFont {Helvetica 12}

	set xstart 100
	set xend [expr $thiswidth-75]
	set ystart [expr $thisheight-75]
	set yend 75

	$c create line $xstart $ystart $xstart $yend -width 2
	$c create line $xstart $ystart $xend $ystart -width 2

	

	
	set unitx [expr double($xend-$xstart)/($xmax-$xmin)]
	set unity [expr double($yend-$ystart)/($ymax-$ymin)]

	for {set i 0} {$i <= 1} {set i [expr $i+0.2]} {
	    $c create line [expr $xstart+$i*($xend-$xstart)] [expr $ystart] [expr $xstart+$i*($xend-$xstart)] [expr $ystart+5] -width 2
	    $c create text [expr $xstart+$i*($xend-$xstart)] [expr $ystart+7] -anchor n -font $smallFont \
		-text [format "%.3g" [expr $xmin+$i*($xmax-$xmin)]]
	    $c create line [expr $xstart] [expr $ystart+$i*($yend-$ystart)] [expr $xstart-7] [expr $ystart+$i*($yend-$ystart)] -width 2
	    $c create text [expr $xstart-9] [expr $ystart+$i*($yend-$ystart)] -anchor e -font $smallFont \
		-text [format "%.3g" [expr $ymin+$i*($ymax-$ymin)]]
	}

	upvar $plotinfo ploti

	set ploti "$xstart $ystart $xmin $ymin $unitx $unity"

		
	button $w.close -text "Close" -command "destroy $w"
	pack $w.close

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w $title
	focus $w

    }

}


proc getlineplotdata { datax datay xmini xmaxi ymini ymaxi} {

    upvar $datax datx
    upvar $datay daty
    
    upvar $xmini xmin
    upvar $xmaxi xmax
    upvar $ymini ymin
    upvar $ymaxi ymax

    global visoptions.lineplotusingx
    global visoptions.lineplotusingy
    global visoptions.lineplotsource
    global visoptions.lineplotfile

    set datx ""
    set daty ""

    set xmin 1e20
    set xmax -1e20
    set ymin 1e20
    set ymax -1e20
    

    if {${visoptions.lineplotsource} == "file"} {
	set fileId [open ${visoptions.lineplotfile} r]
	set line ""
	
	while {[gets $fileId line] >= 0} {
	    if { [string index [lindex $line 0] 0] != "\#" } {
		if { ${visoptions.lineplotusingx} < [llength $line] } {
		    lappend datx [lindex $line ${visoptions.lineplotusingx}]
		    
		    if { [lindex $datx end] < $xmin } {set xmin [lindex $datx end]}
		    if { [lindex $datx end] > $xmax } {set xmax [lindex $datx end]}
		} {
		    lappend datx 0
		}
		if { ${visoptions.lineplotusingy} < [llength $line] } {
		    lappend daty [lindex $line ${visoptions.lineplotusingy}]

		    if { [lindex $daty end] < $ymin } {set ymin [lindex $daty end]}
		    if { [lindex $daty end] > $ymax } {set ymax [lindex $daty end]}
		} {
		    lappend daty 0
		}
	    }
	    
	}
	close $fileId
    }
}


proc lineplotdialog { } {

    set w .lineplot_dlg
    
    if {[winfo exists .lineplot_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	
	toplevel $w
	
	ttk::frame $w.filesettings -relief  groove -borderwidth 3
	ttk::frame $w.filesettings.title
	ttk::radiobutton $w.filesettings.title.choose -variable visoptions.lineplotsource \
	    -value file -text "Data from File"

	pack $w.filesettings.title.choose -side left

	pack $w.filesettings.title

	
	global visoptions.lineplotselectedeval
	global visoptions.lineplotfile
	global visoptions.evaluatefilenames
	global visoptions.evaluatefiledescriptions

	set evdata [NGS_GetData evaluatefiles]
	set visoptions.evaluatefilenames none
	set visoptions.evaluatefiledescriptions none
	for {set i 0} {[expr $i+1] < [llength $evdata]} {incr i 2} {
	    lappend visoptions.evaluatefilenames [lindex $evdata $i]
	    lappend visoptions.evaluatefiledescriptions [lindex $evdata [expr $i+1]]	    
	}
	

	# tixOptionMenu $w.filesettings.latestevals -label "Use Evaluate Results: " \
	    # -options {
		# label.width  25
		# label.anchor e
		# menubutton.width 40
	    # } 
	
	# for {set i 0} {$i < [llength ${visoptions.evaluatefilenames}]} {incr i} {
	    # $w.filesettings.latestevals add command $i \
		# -label "[lindex ${visoptions.evaluatefiledescriptions} $i] ([lindex ${visoptions.evaluatefilenames} $i])"
	# }
	# $w.filesettings.latestevals config -variable visoptions.lineplotselectedeval

	# pack $w.filesettings.latestevals
    ttk::frame $w.filesettings.latestevals
    ttk::label  $w.filesettings.latestevals.lab -text "Use Evaluate Results: "
    ttk::menubutton $w.filesettings.latestevals.but -menu $w.filesettings.latestevals.menu -text "coarse" -width 40

    menu $w.filesettings.latestevals.menu -tearoff 0
	for {set i 0} {$i < [llength ${visoptions.evaluatefilenames}]} {incr i} {
	    $w.filesettings.latestevals.menu add command -label $i\
                -command "set visoptions.lineplotselectedeval $i ; $w.filesettings.latestevals.but configure -text \"[lindex ${visoptions.evaluatefiledescriptions} $i] ([lindex ${visoptions.evaluatefilenames} $i])\""
	}
   $w.filesettings.latestevals.menu invoke ${visoptions.lineplotselectedeval}               


    grid $w.filesettings.latestevals.lab $w.filesettings.latestevals.but -sticky nw
	pack $w.filesettings.latestevals
	ttk::frame $w.filesettings.sfn

	ttk::button $w.filesettings.sfn.bb -text "Browse" \
	    -command { set visoptions.lineplotfile [tk_getOpenFile] }

	
	ttk::entry $w.filesettings.sfn.fn -width 50 \
	    -textvariable visoptions.lineplotfile

	pack $w.filesettings.sfn.bb $w.filesettings.sfn.fn -side left

	pack $w.filesettings.sfn

	ttk::button $w.filesettings.refresh -text "Refresh" -command {
	    if { ${visoptions.lineplotselectedeval} != 0} {
		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]
	    }
	    
	    set saveusingx ${visoptions.lineplotusingx}
	    set saveusingy ${visoptions.lineplotusingy}
	    

	    
	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
		${visoptions.lineplotxcoordselector} delete $i
	    }
	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
		${visoptions.lineplotycoordselector} delete $i
	    }

	    
	    set fileId [open ${visoptions.lineplotfile} r]
	    set line ""
	    gets $fileId line
	    close $fileId
	    if { [lindex $line 0] == "\#nglineplotinfo" } {
		set visoptions.lineplotdatadescr [lrange $line 1 end]	
	    } {
		set visoptions.lineplotdatadescr ""
		for { set i 0 } { $i < [llength $line] } { incr i } {
		    lappend visoptions.lineplotdatadescr "data[expr $i+1]"
		}
	    }

	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
		${visoptions.lineplotxcoordselector} add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	    }
	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
		${visoptions.lineplotycoordselector} add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	    }

	    if { $saveusingx < [llength ${visoptions.lineplotdatadescr}] } {
		set visoptions.lineplotusingx $saveusingx
	    } {
		set visoptions.lineplotusingx 0
	    }
	    if { $saveusingy < [llength ${visoptions.lineplotdatadescr}] } {
		set visoptions.lineplotusingy $saveusingy
	    } {
		set visoptions.lineplotusingy 1
	    }	    
	}

	pack $w.filesettings.refresh


	ttk::frame $w.filesettings.using

	global visoptions.lineplotdatadescr

	# tixOptionMenu $w.filesettings.using.xco -label "X-Coord:"\
	    # -options {
		# label.width  8
		# label.anchor e
		# menubutton.width 15
	    # } 
	# for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
	    # $w.filesettings.using.xco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	# }
    
    ttk::frame $w.filesettings.using.xco
    ttk::label  $w.filesettings.using.xco.lab -text "X-Coord:"
    ttk::menubutton $w.filesettings.using.xco.but -menu $w.filesettings.using.xco.menu -text "" -width 15

    menu $w.filesettings.using.xco.menu -tearoff 0
	for {set i 0} {$i < [llength ${visoptions.lineplotdatadescr}]} {incr i} {
	    $w.filesettings.using.xco.menu add command -label [lindex ${visoptions.lineplotdatadescr} $i]\
                -command "set visoptions.lineplotusingx $i ; $w.filesettings.using.xco.but configure -text \"[lindex ${visoptions.lineplotdatadescr} $i]\""
	}
   $w.filesettings.using.xco.menu invoke [lindex ${visoptions.lineplotdatadescr} 0]


    grid $w.filesettings.using.xco.lab $w.filesettings.using.xco.but -sticky nw
	#pack $w.filesettings.using.xco    
    
    
    
	# $w.filesettings.using.xco config -variable visoptions.lineplotusingx
	
	# tixOptionMenu $w.filesettings.using.yco -label "Y-Coord:"\
	    # -options {
		# label.width  8
		# label.anchor e
		# menubutton.width 15
	    # } 
	# for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
	    # $w.filesettings.using.yco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	# }
	# $w.filesettings.using.yco config -variable visoptions.lineplotusingy
    ttk::frame $w.filesettings.using.yco
    ttk::label  $w.filesettings.using.yco.lab -text "Y-Coord:"
    ttk::menubutton $w.filesettings.using.yco.but -menu $w.filesettings.using.yco.menu -text "" -width 15

    menu $w.filesettings.using.yco.menu -tearoff 0
	for {set i 0} {$i < [llength ${visoptions.lineplotdatadescr}]} {incr i} {
	    $w.filesettings.using.yco.menu add command -label [lindex ${visoptions.lineplotdatadescr} $i]\
                -command "set visoptions.lineplotusingy $i ; $w.filesettings.using.yco.but configure -text \"[lindex ${visoptions.lineplotdatadescr} $i]\""
	}
   $w.filesettings.using.yco.menu invoke [lindex ${visoptions.lineplotdatadescr} 0]
   grid $w.filesettings.using.yco.lab $w.filesettings.using.yco.but -sticky nw
    
	global visoptions.lineplotxcoordselector
	global visoptions.lineplotycoordselector
	set visoptions.lineplotxcoordselector $w.filesettings.using.xco
	set visoptions.lineplotycoordselector $w.filesettings.using.yco
	

	pack $w.filesettings.using.xco $w.filesettings.using.yco -side left
	pack $w.filesettings.using

	pack $w.filesettings -fill x -ipady 3
	
	ttk::frame $w.settings -relief  groove -borderwidth 3
	ttk::label $w.settings.title -text "\nSettings\n"
	pack $w.settings.title 

	ttk::frame $w.settings.minmax 
	ttk::checkbutton $w.settings.minmax.autoscale -text "Autoscale" -variable visoptions.lineplotautoscale
	# tixControl $w.settings.minmax.xmin -label "Min. x: " \
	    # -integer false -variable visoptions.lineplotxmin \
	    # -options {
		# entry.width 6
		# label.width 8
		# label.anchor e
	    # }	
    ttk::frame $w.settings.minmax.xmin
    ttk::label $w.settings.minmax.xmin.label -text "Min. x: "
    ttk::spinbox $w.settings.minmax.xmin.sp -textvariable visoptions.lineplotxmin -width 6 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 3" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9 

        # tixControl $w.settings.minmax.xmax -label "Max. x: " \
	    # -integer false -variable visoptions.lineplotxmax \
	    # -options {
		# entry.width 6
		# label.width 8
		# label.anchor e
	    # }	
    ttk::frame $w.settings.minmax.xmax
    ttk::label $w.settings.minmax.xmax.label -text "Max. x: "
    ttk::spinbox $w.settings.minmax.xmax.sp -textvariable visoptions.lineplotxmax -width 6 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 3" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9 
        
	# tixControl $w.settings.minmax.ymin -label "Min. y: " \
	    # -integer false -variable visoptions.lineplotymin \
	    # -options {
		# entry.width 6
		# label.width 8
		# label.anchor e
	    # }	
    ttk::frame $w.settings.minmax.ymin
    ttk::label $w.settings.minmax.ymin.label -text "Min. y: "
    ttk::spinbox $w.settings.minmax.ymin.sp -textvariable visoptions.lineplotymin -width 6 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 3" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9 
        
	# tixControl $w.settings.minmax.ymax -label "Max. y: " \
	    # -integer false -variable visoptions.lineplotymax \
	    # -options {
		# entry.width 6
		# label.width 8
		# label.anchor e
	    # }	
    ttk::frame $w.settings.minmax.ymax
    ttk::label $w.settings.minmax.ymax.label -text "Max. y: "
    ttk::spinbox $w.settings.minmax.ymax.sp -textvariable visoptions.lineplotymax -width 6 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 3" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9 
    pack $w.settings.minmax.xmin.label $w.settings.minmax.xmin.sp
    pack $w.settings.minmax.xmax.label $w.settings.minmax.xmax.sp
    pack $w.settings.minmax.ymin.label $w.settings.minmax.ymin.sp
    pack $w.settings.minmax.ymax.label $w.settings.minmax.ymax.sp
	pack $w.settings.minmax.autoscale $w.settings.minmax.xmin $w.settings.minmax.xmax \
	    $w.settings.minmax.ymin $w.settings.minmax.ymax -side left

	pack $w.settings.minmax

	
	ttk::label $w.settings.empty1 -text ""
	pack $w.settings.empty1

	ttk::frame $w.settings.plotsize

	# tixControl $w.settings.plotsize.xsize -label "Plotsize  x: "\
	    # -integer true -variable visoptions.lineplotsizex \
	    # -options {
		# entry.width 6
		# label.width 13
		# label.anchor e
	    # }
        
    ttk::frame $w.settings.plotsize.xsize
    ttk::label $w.settings.plotsize.xsize.label -text "Plotsize  x: "
    ttk::spinbox $w.settings.plotsize.xsize.sp -textvariable visoptions.lineplotsizex -width 6 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9
    pack $w.settings.plotsize.xsize.label $w.settings.plotsize.xsize.sp
    
	# tixControl $w.settings.plotsize.ysize -label "y: "\
	    # -integer true -variable visoptions.lineplotsizey \
	    # -options {
		# entry.width 6
		# label.width 3
		# label.anchor e
	    # }	
    
    ttk::frame $w.settings.plotsize.ysize
    ttk::label $w.settings.plotsize.ysize.label -text "Plotsize  y: "
    ttk::spinbox $w.settings.plotsize.ysize.sp -textvariable visoptions.lineplotsizey -width 6 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
        -invalidcommand "my_invalidspinbox %W" -from -1e9 -to 1e9
    pack $w.settings.plotsize.ysize.label $w.settings.plotsize.ysize.sp

	pack $w.settings.plotsize.xsize $w.settings.plotsize.ysize -side left

	pack $w.settings.plotsize

	ttk::label $w.settings.empty2 -text ""
	pack $w.settings.empty2
	
	# tixOptionMenu $w.settings.color -label "Linecolor: " \
	    # -options {
		# label.width  19
		# label.anchor e
		# menubutton.width 15
	    # }
	# foreach step { red black blue green yellow } {
	    # $w.settings.color add command $step -label $step
	# }
	# $w.settings.color config -variable visoptions.lineplotcolor
    ttk::frame $w.settings.color
    ttk::label  $w.settings.color.lab -text "Linecolor: "
    ttk::menubutton $w.settings.color.but -menu $w.settings.color.menu -text "" -width 15

    menu $w.settings.color.menu -tearoff 0
	foreach step { red black blue green yellow } {
	    $w.settings.color.menu add command -label $step -command "set visoptions.lineplotcolor $step; $w.settings.color.but configure -text \"$step\""
    }
	# for {set i 0} {$i < [llength ${visoptions.lineplotdatadescr}]} {incr i} {
	    # $w.filesettings.using.yco.menu add command -label [lindex ${visoptions.lineplotdatadescr} $i]\
                # -command "set visoptions.lineplotusingy $i ; $w.filesettings.using.yco.but configure -text \"[lindex ${visoptions.lineplotdatadescr} $i]\""
	# }
   $w.settings.color.menu invoke "red"
   grid $w.settings.color.lab $w.settings.color.but -sticky nw    
    
	pack $w.settings.color


	pack $w.settings -fill x 

	set datax ""
	set datay ""
	set xmin 0
	set xmax 0
	set ymin 0
	set ymax 0

	ttk::frame $w.plots -relief  groove -borderwidth 3

	# tixOptionMenu $w.plots.selplot -label "Selected Plot: " \
	    # -options {
		# label.width  19
		# label.anchor e
		# menubutton.width 15
	    # } 
	# $w.plots.selplot add command none -label "None"

	# $w.plots.selplot config -variable visoptions.lineplotselected

    ttk::frame $w.plots.selplot
    ttk::label  $w.plots.selplot.lab -text "Linecolor: "
    ttk::menubutton $w.plots.selplot.but -menu $w.plots.selplot.menu -text "" -width 15

    menu $w.plots.selplot.menu -tearoff 0
	$w.plots.selplot.menu add command -label "None" -command "set visoptions.lineplotselected \"None\"; $w.plots.selplot.but configure -text \"None\""
    grid $w.plots.selplot.lab $w.plots.selplot.but -sticky nw    
    $w.plots.selplot.menu invoke "None"
    
	global visoptions.lineplotselector
	set visoptions.lineplotselector $w.plots.selplot.menu

	

	ttk::button $w.plots.new -text "Generate New Plot" -command {
	    if { ${visoptions.lineplotselectedeval} != 0} {
		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]
	    }

	    getlineplotdata datax datay xmin xmax ymin ymax

	    puts stdout "xmin $xmin xmax $xmax ymin $ymin ymax $ymax"

	    global visoptions.lineplotautoscale

	    if {! ${visoptions.lineplotautoscale}} {
		puts "using set min/max values"
		set xmin ${visoptions.lineplotxmin}
		set xmax ${visoptions.lineplotxmax}
		set ymin ${visoptions.lineplotymin}
		set ymax ${visoptions.lineplotymax}
	    }
	    
	    incr visoptions.lineplotcurrentnum
	    
	    set ident .newplot${visoptions.lineplotcurrentnum}
	    set plotinfo ""
	    
	    createlineplot ${visoptions.lineplotsizex} ${visoptions.lineplotsizey} \
		$ident "Lineplot ${visoptions.lineplotcurrentnum}" \
		$xmin $xmax $ymin $ymax plotinfo
	    
	    lappend visoptions.lineplotinfos $plotinfo

	    
	    ${visoptions.lineplotselector} add command ${visoptions.lineplotcurrentnum} -label "Lineplot ${visoptions.lineplotcurrentnum}"

	    addplotline $ident $datax $datay $plotinfo ${visoptions.lineplotcolor}
	}
	
	ttk::button $w.plots.addto -text "Add to Selected Plot" -command {
	    if { ${visoptions.lineplotselectedeval} != 0} {
		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]
	    }

	    if { ${visoptions.lineplotselected} != "none" } {

		getlineplotdata datax datay xmin xmax ymin ymax

		set ident .newplot${visoptions.lineplotselected}
		set plotinfo [lindex ${visoptions.lineplotinfos} ${visoptions.lineplotselected}]

		addplotline $ident $datax $datay $plotinfo ${visoptions.lineplotcolor}
	    }
	}
	    
	
	
	pack $w.plots.new $w.plots.addto $w.plots.selplot
	
	
	pack $w.plots -fill x -ipady 3
	


	ttk::button $w.close -text "Close" -command "destroy $w"
	pack $w.close
	
	wm withdraw $w
	wm geom $w +200+100
	wm deiconify $w
	wm title $w "2D Lineplots"

	focus $w
    

    }
}


set fieldlinesdialog_pop1 0


proc fieldlinesdialog { } {
    
    set w .fieldlines_dlg

    global fieldlinesdialog_pop1
    set fieldlinesdialog_pop1 1

    
    if {[winfo exists .fieldlines_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	
	toplevel $w

	#tixNoteBook $w.nb -ipadx 6 -ipady 6

        pack [ttk::notebook $w.nb]  -fill both -side top -ipadx 6 -ipady 6
        $w.nb add [ttk::frame $w.nb.draw] -text "Draw" -underline 0 
        $w.nb add [ttk::frame $w.nb.settings] -text "Settings" -underline 0        
        
	#$w.nb add draw -label "Draw"
	#$w.nb add settings -label "Settings"

	#pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top


	# Main Window

	set f $w.nb.draw
	
	ttk::labelframe $f.general -text "General settings" -relief groove -borderwidth 3

	
	ttk::checkbutton $f.general.enable -text "Enable Fieldlines" \
	    -variable visoptions.drawfieldlines \
	    -command { 
		# set visoptions.redrawperiodic ${visoptions.drawfieldlines}
		# redrawperiodic
		# redrawperiodic # sonst 
		Ng_Vis_Set parameters; 
		redraw 
	    }
        ttk::label $f.general.numl -text "num:"
	ttk::spinbox $f.general.num -from 0 -to 100 -increment 1 \
            -textvariable visoptions.numfieldlines -width 4
        #tixControl $f.general.num -label "Num: " -integer true \
	    -variable visoptions.numfieldlines \
	    -command { Ng_Vis_Set parameters; redraw } \
	    -options {
	#	entry.width 6
	#	label.width 12
	#	label.anchor e
	#    }	

	
	grid $f.general.enable -sticky nw -padx 4 -pady 2
        grid x $f.general.numl $f.general.num -rowspan 3 -sticky w -padx 4 -row 0 -pady 2
        grid anchor $f.general center

	pack $f.general -pady 15 -fill x -ipady 3

	ttk::label $f.labe0 -text " "

	pack $f.labe0

	#ttk::frame $f.general1
	
	ttk::checkbutton $f.general.randomstart -text "Field dependent density    " \
	    -variable visoptions.fieldlinesrandomstart \
	    -command { Ng_Vis_Set parameters; redraw}

	
	ttk::checkbutton $f.general.redrawperiodic -text "Animate periodic" \
	    -variable visoptions.redrawperiodic \
	    -command { 
		redrawperiodic
		Ng_Vis_Set parameters; 
		redraw 
	    }

	grid $f.general.randomstart -sticky nw -padx 4 -row 1
        grid $f.general.redrawperiodic -sticky nw -padx 4 -row 2

	#pack $f.general1



	ttk::label $f.lab0 -text " "

	pack $f.lab0


	
	# tixOptionMenu $f.vecfun -label "Vector Function: " \
	    # -options {
		# label.width  18
		# label.anchor e
		# menubutton.width 12
	    # }
	# $f.vecfun add command none -label None 
	# for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    # set fname [Ng_Vis_Field getfieldname $i]
	    # set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    # set iscomplex [Ng_Vis_Field iscomplex $i]
	    # set sdim [Ng_Vis_Field getdimension]
	    # if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    # if { ($fcomp == $sdim) || ($fcomp == 3) } {
		# $f.vecfun add command $fname -label $fname
	    # } 
	# }
	# $f.vecfun configure -variable visoptions.fieldlinesvecfunction
	# $f.vecfun configure -command { Ng_Vis_Set parameters; redraw }
    ttk::frame $f.vecfun
    ttk::label  $f.vecfun.lab -text "Vector Function: "     
    ttk::menubutton $f.vecfun.but -menu $f.vecfun.menu -text "" -width 12     
    menu $f.vecfun.menu -tearoff 0 	
    # for {set i 0} {$i < [llength ${visoptions.evaluatefilenames}]} {incr i} {
    # $w.filesettings.latestevals.menu add command -label $i\                 
    # -command "set visoptions.lineplotselectedeval $i ; $w.filesettings.latestevals.but configure -text \"[lindex ${visoptions.evaluatefiledescriptions} $i] ([lindex ${visoptions.evaluatefilenames} $i])\"" 	
    # } 	
    for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
    set fname [Ng_Vis_Field getfieldname $i] 	    
    set fcomp [Ng_Vis_Field getfieldcomponents $i] 	    
    set iscomplex [Ng_Vis_Field iscomplex $i] 	    
    set sdim [Ng_Vis_Field getdimension] 	    
    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] } 	    
    if { ($fcomp == $sdim) || ($fcomp == 3) } {
    $f.vecfun.menu add command -label $fname -command "set visoptions.fieldlinesvecfunction $fname;Ng_Vis_Set parameters; redraw;$f.vecfun.but configure -text \"$fname\" " 	    }         
    }
    grid $f.vecfun.lab $f.vecfun.but -sticky nw
	pack $f.vecfun
	


	ttk::label $f.lab00 -text " "

	pack $f.lab00

	ttk::frame $f.phasesettings

	ttk::checkbutton $f.phasesettings.onephase -text "Fix Phase" -variable visoptions.fieldlinesonlyonephase
	# scale $f.phasesettings.phase -orient horizontal -length 300 -from 0 -to 360 \
	    # -label "phi" \
	    # -resolution 1 \
	    # -variable visoptions.fieldlinesphase \
	    # -command { popupcheckredraw3 fieldlinesdialog_pop1 }

    ttk::frame $f.phasesettings.phase
    ttk::label $f.phasesettings.phase.lab -text "phi"
    ttk::scale $f.phasesettings.phase.sc -orient horizontal -length 100 -from 0 -to 360 -variable visoptions.fieldlinesphase \
        -command "roundscale $f.phasesettings.phase.sc 0; popupcheckredraw3 fieldlinesdialog_pop1"
    ttk::entry $f.phasesettings.phase.ent -width 4 -textvariable visoptions.fieldlinesphase -validate focus -takefocus 0 \
        -validatecommand "popupcheckredraw3 fieldlinesdialog_pop1;my_validate %W 0 1 %P 0" \
        -invalidcommand "my_invalid %W;popupcheckredraw3 fieldlinesdialog_pop1"
    grid $f.phasesettings.phase.lab $f.phasesettings.phase.sc $f.phasesettings.phase.ent -sticky nw -ipadx 4
	pack $f.phasesettings.onephase $f.phasesettings.phase -side left
	
	pack $f.phasesettings



	ttk::label $f.lab1 -text " "

	pack $f.lab1


	
	ttk::frame $f.boxsettings -relief groove -borderwidth 3
	ttk::frame $f.boxsettings.title
	ttk::radiobutton $f.boxsettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value box -text "Startpoints in Box"

	pack $f.boxsettings.title.choose -side left

	pack $f.boxsettings.title

	ttk::frame $f.boxsettings.points

	ttk::label $f.boxsettings.points.lab2 -text "Pmin";
	ttk::entry $f.boxsettings.points.ent1x -width 8 \
	    -textvariable visoptions.fieldlinesstartareap1x
	ttk::entry $f.boxsettings.points.ent1y -width 8 \
	    -textvariable visoptions.fieldlinesstartareap1y
	ttk::entry $f.boxsettings.points.ent1z -width 8 \
	    -textvariable visoptions.fieldlinesstartareap1z
	ttk::label $f.boxsettings.points.lab3 -text "   Pmax";
	ttk::entry $f.boxsettings.points.ent2x -width 8 \
	    -textvariable visoptions.fieldlinesstartareap2x
	ttk::entry $f.boxsettings.points.ent2y -width 8 \
	    -textvariable visoptions.fieldlinesstartareap2y
	ttk::entry $f.boxsettings.points.ent2z -width 8 \
	    -textvariable visoptions.fieldlinesstartareap2z
	
	pack $f.boxsettings.points
	pack $f.boxsettings.points.lab2 $f.boxsettings.points.ent1x $f.boxsettings.points.ent1y $f.boxsettings.points.ent1z -side left
	pack $f.boxsettings.points.lab3 $f.boxsettings.points.ent2x $f.boxsettings.points.ent2y $f.boxsettings.points.ent2z -side left

	ttk::button $f.boxsettings.settobb -text "Bounding Box" -command {
	    set bbox [Ng_MeshInfo bbox]
	    set visoptions.fieldlinesstartareap1x [lindex $bbox 0]
	    set visoptions.fieldlinesstartareap2x [lindex $bbox 1]
	    set visoptions.fieldlinesstartareap1y [lindex $bbox 2]
	    set visoptions.fieldlinesstartareap2y [lindex $bbox 3]
	    set visoptions.fieldlinesstartareap1z [lindex $bbox 4]
	    set visoptions.fieldlinesstartareap2z [lindex $bbox 5]
	}

	pack $f.boxsettings.settobb
	    
	
	pack $f.boxsettings -fill x -ipady 3


	ttk::frame $f.facesettings -relief groove -borderwidth 3
	ttk::frame $f.facesettings.title
	ttk::radiobutton $f.facesettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value face -text "Startpoints on Face"

	pack $f.facesettings.title.choose -side left

	pack $f.facesettings.title
	
	ttk::frame $f.facesettings.index
	ttk::label $f.facesettings.index.lab -text "face index:"
	ttk::label $f.facesettings.index.ent -text 1;# -padx 4

	pack $f.facesettings.index.lab $f.facesettings.index.ent -side left

	pack $f.facesettings.index
	
	pack $f.facesettings -fill x -ipady 3


	global visoptions.fieldlinesfilename

	ttk::frame $f.filesettings -relief  groove -borderwidth 3
	ttk::frame $f.filesettings.title
	ttk::radiobutton $f.filesettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value file -text "Startpoints from File"

	pack $f.filesettings.title.choose -side left

	pack $f.filesettings.title

	ttk::frame $f.filesettings.sfn

	ttk::button $f.filesettings.sfn.bb -text "Browse" \
	    -command {
		set types {
		    { "Netgen Fieldlines" {.nef} }
		}
		set visoptions.fieldlinesfilename [tk_getOpenFile -filetypes $types -defaultextension ".nef"]
	    }

	
	ttk::entry $f.filesettings.sfn.fn -width 50 \
	    -textvariable visoptions.fieldlinesfilename

	pack $f.filesettings.sfn.bb $f.filesettings.sfn.fn -side left

	pack $f.filesettings.sfn

	pack $f.filesettings -fill x -ipady 3
	


	
	# Settings

	set g $w.nb.settings

	ttk::frame $g.linesettings -relief groove -borderwidth 3
	ttk::label $g.linesettings.title -text "\nLine Settings\n"
	# tixControl $g.linesettings.length -label "rel. Length: " -integer false \
	    # -variable visoptions.fieldlineslength -min 0.00001 -max 10000 -step 0.1 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }

    ttk::frame $g.linesettings.length
    ttk::label $g.linesettings.length.lab -text "rel. Length: "
    ttk::spinbox $g.linesettings.length.sp -textvariable visoptions.fieldlineslength -width 6 -increment 0.1 -validate focus -validatecommand "my_validatespinbox %W %P 5" \
        -invalidcommand "my_invalidspinbox %W" -from 0.00001 -to 10000
    grid $g.linesettings.length.lab $g.linesettings.length.sp -sticky nw
        
	# tixControl $g.linesettings.maxpoints -label "max. Points: " -integer true \
	    # -variable visoptions.fieldlinesmaxpoints -min 0 -max 10000 -step 1 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }
    ttk::frame $g.linesettings.maxpoints
    ttk::label $g.linesettings.maxpoints.lab -text "max. Points: "
    ttk::spinbox $g.linesettings.maxpoints.sp -textvariable visoptions.fieldlinesmaxpoints -width 6 -increment 1 -validate focus -validatecommand "my_validatespinbox %W %P 0" \
        -invalidcommand "my_invalidspinbox %W" -from 0 -to 10000
    grid $g.linesettings.maxpoints.lab $g.linesettings.maxpoints.sp -sticky nw

        
	# tixControl $g.linesettings.thick -label "rel. Thickness: " -integer false \
	    # -variable visoptions.fieldlinesthickness -min 1e-10 -max 0.5 -step 0.001 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }
        
    ttk::frame $g.linesettings.thick
    ttk::label $g.linesettings.thick.lab -text "rel. Thickness: "
    ttk::spinbox $g.linesettings.thick.sp -textvariable visoptions.fieldlinesthickness -width 6 -increment 0.001 -validate focus -validatecommand "my_validatespinbox %W %P 6" \
        -invalidcommand "my_invalidspinbox %W" -from 1e-10 -to 0.5
    grid $g.linesettings.thick.lab $g.linesettings.thick.sp -stick nw
    
	pack $g.linesettings.title $g.linesettings.length $g.linesettings.maxpoints $g.linesettings.thick

	pack $g.linesettings -fill x -ipady 3


	


	global visoptions.fieldlinestolerance

	ttk::frame $g.odesettings -relief groove -borderwidth 3
	ttk::label $g.odesettings.title -text "\nODE Settings\n"
	# tixControl $g.odesettings.tol -label "rel. Tolerance: " -integer false \
	    # -variable visoptions.fieldlinestolerance -min 0.00001 -max 1 -step 0.01 \
	    # -options {
		# entry.width 6
		# label.width 25
		# label.anchor e
	    # }	
    
    ttk::frame $g.odesettings.tol
    ttk::label $g.odesettings.tol.lab -text "rel. Thickness: "
    ttk::spinbox $g.odesettings.tol.sp -textvariable visoptions.fieldlinestolerance -width 6 -increment 0.01 -validate focus -validatecommand "my_validatespinbox %W %P 5" \
        -invalidcommand "my_invalidspinbox %W" -from 0.00001 -to 1
    grid $g.odesettings.tol.lab $g.odesettings.tol.sp -stick nw

        
	# tixOptionMenu $g.odesettings.rktype -label "RK-Type " \
	    # -options {
		# label.width  20
		# label.anchor e
		# menubutton.width 25
	    # }
	# $g.odesettings.rktype add command euler -label "Euler, order 1"
	# $g.odesettings.rktype add command eulercauchy -label "Euler-Cauchy, order 2"
	# $g.odesettings.rktype add command simpson -label "Simpson, order 3"
	# $g.odesettings.rktype add command crungekutta -label "classical Runge-Kutta, order 4"
	# $g.odesettings.rktype configure -variable visoptions.fieldlinesrktype
	# $g.odesettings.rktype configure -command { Ng_Vis_Set parameters; redraw }
    ttk::frame $g.odesettings.rktype     
    ttk::label  $g.odesettings.rktype.lab -text "RK-Type "     
    ttk::menubutton $g.odesettings.rktype.but -menu $g.odesettings.rktype.menu -text "" -width 25     
    menu $g.odesettings.rktype.menu -tearoff 0     
    $g.odesettings.rktype.menu add command -label "Euler, order 1" -command "set visoptions.fieldlinesrktype \"euler\" ;Ng_Vis_Set parameters; redraw;$g.odesettings.rktype.but configure -text \"Euler,order 1\" "     
    $g.odesettings.rktype.menu add command -label "Euler-Cauchy, order 2" -command "set visoptions.fieldlinesrktype \"eulercauchy\" ;Ng_Vis_Set parameters; redraw;$g.odesettings.rktype.but configure -text \"Euler-Cauchy,order 2\" "     
    $g.odesettings.rktype.menu add command -label "Simpson, order 3" -command "set visoptions.fieldlinesrktype \"simpson\" ;Ng_Vis_Set parameters; redraw;$g.odesettings.rktype.but configure -text \"Simpson,order 3\""     
    $g.odesettings.rktype.menu add command -label "classical Runge-Kutta, order 4" -command "set visoptions.fieldlinesrktype \"crungekutta\" ;Ng_Vis_Set parameters; redraw; $g.odesettings.rktype.but configure -text \"classical Runge-Kutta,order 4\""        
    $g.odesettings.rktype.menu invoke "classical Runge-Kutta, order 4"     
    
    grid $g.odesettings.rktype.lab $g.odesettings.rktype.but -sticky nw 	
	
	pack $g.odesettings.title $g.odesettings.tol $g.odesettings.rktype

	pack $g.odesettings -fill x -ipady 3



	# buttons
	

	ttk::frame $w.bu 
	pack $w.bu -fill x -ipady 3

	ttk::button $w.bu.calc -text "Build Fieldlines" -command { 
	    if { ${visoptions.fieldlinesvecfunction} == "none" } {
		bgerror "Please select the vector function first!"
	    } {
		set visoptions.drawfieldlines 1
		Ng_Vis_Set parameters
		Ng_BuildFieldLines
		redraw 
	    }
	}

	ttk::button $w.bu.help -text "Help" -command {
	    if {[winfo exists .fieldlines_help] == 1} {
		wm withdraw .fieldlines_help
		wm deiconify .fieldlines_help
		focus .fieldlines_help
	    } {
		toplevel .fieldlines_help

        set f [frame .fieldlines_help.ht]
        #ttk::scrollbar $f.hsb -orient horizontal -command [list $f.t xview]
        ttk::scrollbar $f.vsb -orient vertical -command [list $f.t yview]
        text $f.t -yscrollcommand [list $f.vsb set] 
        grid $f.t -row 0 -column 0 -sticky nsew
        grid $f.vsb -row 0 -column 1 -sticky nsew
        grid columnconfigure $f 0 -weight 1
        grid rowconfigure $f 0 -weight 1
		#tixScrolledText .fieldlines_help.ht -scrollbar y
		set text $f.t

		$text configure -setgrid true -wrap word 

		$text tag configure bold -font *-*-bold-*-*-*-*

		
		$text insert end \
		    "Draw menu\n \n" bold
		$text insert end \
		    "Enable Fieldlines\n    To turn on and off the calculated fieldlines. (Has to be turned on to start the calculation)\n"
		$text insert end \
		    "Num\n    Number of fieldlines to calculate. (May not be used exactly.)"
		$text insert end \
		    "Field dependent density\n    There will be more fieldline startpoints where the field is stronger\n\n"
		$text insert end \
		    "Animate periodic\n    (for quasistationary fields) The fieldlines of the different phase angles are animated.\n    ATTENTION: \"Fix Phase\" has to be turned off\n\n"
		$text insert end \
		    "Vector Function\n    The function fixing the direction of the lines\n\n"
		$text insert end \
		    "Fix Phase\n    (for quasistationary fields) Only calculate and draw fieldlines for one special phase angle.\n\n"
		$text insert end \
		    "Startpoints in Box\n    Set the startpoints inside the box \[Pmin1,Pmax1\] x \[Pmin2,Pmax2\] x \[Pmin3,Pmax3\]\n"
		$text insert end \
		    "    With the button \"Bounding Box\" the whole bounding box of the geometry is selected.\n\n" 
		$text insert end \
		    "Startpoints on Face\n    All startpoints will be set on one face. This face is selected by double-clicking with the mouse.\n\n"
		$text insert end \
		    "Startpoints from File\n    The startpoint information will be read from the selected file.\n    The entries in the file can be as follows:\n"
		$text insert end \
		    "        point <x> <y> <z>\n            set a (potential) startpoint\n"
		$text insert end \
		    "        line <x1> <y1> <z1> <x2> <y2> <z2> <n>\n            set n (potential) startpoints on the line from (x1,y1,z1) to (x2,y2,z2)\n"
		$text insert end \
		    "        box <x1> <y1> <z1> <x2> <y2> <z2> <n>\n            set n (potential) startpoints inside the box \[x1,x2\] x \[y1,y2\] x \[z1,z2\]\n"
		$text insert end \
		    "    ATTENTION: These are potential startpoints.\n               The total number of startpoints will be bounded by the \"Num\"-parameter.\n \n \n \n"
		$text insert end \
		    "Settings Menu\n \n" bold
		$text insert end \
		    "rel. Length\n    The maximal length of a fieldline relative to the diameter of the geometry.\n\n"
		$text insert end \
		    "max. Points\n    The maximum number of Runge-Kutta steps.\n\n"
		$text insert end \
		    "rel. Thickness\n    The thickness of the fieldlines relative to the diameter of the geometry.\n\n"
		$text insert end \
		    "rel. Tolerance\n    The tolerance for the step-length control of the Runge-Kutta method.\n\n"
		$text insert end \
		    "RK-Type\n    Which Runge-Kutta scheme to use\n \n \n \n"
		$text insert end \
		    "Button \"Build Fieldlines\"\n" bold
		$text insert end \
		    "    Build the fieldlines."
		

		$text configure -state disabled

		pack .fieldlines_help.ht -expand yes -fill both

		wm withdraw .fieldlines_help
		wm geom .fieldlines_help +300+200
		wm deiconify .fieldlines_help
		wm title .fieldlines_help "Fieldlines Help"
		focus .fieldlines_help
		
	    }


	}

	ttk::button $w.bu.cancel -text "Done" -command "destroy $w"
	grid $w.bu.calc $w.bu.help $w.bu.cancel -sticky nw -padx 4
        grid anchor $w.bu center
	
	
	wm withdraw $w
	wm geom $w +200+100
	wm deiconify $w
	wm title $w "Fieldlines"
	#    grab $w
	focus $w

    }

    global visoptions.fieldlinesstartface

    
    set f $w.nb.draw
    set visoptions.fieldlinesstartface [Ng_BCProp getactive]
    $f.facesettings.index.ent configure -text ${visoptions.fieldlinesstartface}

}


#proc popupcheckredraw { vari { x 0 } } {
#    upvar $vari varname
#    if { $varname == 1 } {
#	set varname 0
#    } {
#	Ng_Vis_Set parameters
#	redraw
#    }
#}


set visual_dialog_pop1 0
set visual_dialog_pop2 0
set visual_dialog_pop3 0
set visual_dialog_pop4 0
set visual_dialog_pop5 0
set visual_dialog_pop6 0
set visual_dialog_pop7 0

proc visual_dialog { } {

    set w .visoptions_dlg

    global visual_dialog_pop1
    global visual_dialog_pop2
    global visual_dialog_pop3
    global visual_dialog_pop4
    global visual_dialog_pop5
    global visual_dialog_pop6
    global visual_dialog_pop7
    set visual_dialog_pop1 1
    set visual_dialog_pop2 1
    set visual_dialog_pop3 1
    set visual_dialog_pop4 1
    set visual_dialog_pop5 1
    set visual_dialog_pop6 1
    set visual_dialog_pop7 1
   
    if {[winfo exists .visoptions_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w
    } {

	toplevel $w 
        
	#ttk::frame $w.grid -relief groove -borderwidth 0
	# change to: max gridsize 200
	#scale $w.grid.size -orient horizontal -length 100 -from 1 -to 200 \
	    -label "Grid" \
	    -resolution 1    \
	    -variable  visoptions.gridsize \
	    -command { popupcheckredraw visual_dialog_pop2 }
        

	# x- and y- offset
	#scale $w.grid.xoffset -orient horizontal -length 80 -from 0 -to 1 \
	    -label "x-Offset" \
	    -resolution 0.05    \
	    -variable  visoptions.xoffset \
	    -command { popupcheckredraw visual_dialog_pop3 }
        ttk::frame $w.main
        pack $w.main -fill x 
        set w $w.main
        ttk::frame $w.upperfr ;# -relief groove -borderwidth 3 -height 10
        pack $w.upperfr -fill x;# -ipady 8
        
        ttk::labelframe $w.upperfr.size -text "Grid" -relief groove -borderwidth 3
        ttk::entry $w.upperfr.size.ent -width 3 -textvariable visoptions.gridsize -validate focus -takefocus 0 -validatecommand "popupcheckredraw visual_dialog_pop2;my_validate %W 0 200 %P 0" \
            -invalidcommand "my_invalid %W;popupcheckredraw visual_dialog_pop2"
        ttk::scale $w.upperfr.size.sc -orient horizontal -length 100 -from 1 -to 200 -variable visoptions.gridsize\
            -command "roundscale $w.upperfr.size.sc 0;popupcheckredraw visual_dialog_pop2"

        ttk::labelframe $w.upperfr.offsets -text "x / y offsets" -relief groove -borderwidth 3
        ttk::label $w.upperfr.offsets.xlab -text "x"
        ttk::label $w.upperfr.offsets.ylab -text "y"
        ttk::scale $w.upperfr.offsets.xoffset -orient horizontal -length 100 -from 0 -to 1 -variable visoptions.xoffset \
            -command "roundscale $w.upperfr.offsets.xoffset 2; popupcheckredraw visual_dialog_pop3"
        ttk::scale $w.upperfr.offsets.yoffset -orient horizontal -length 100 -from 0 -to 1 -variable visoptions.yoffset \
            -command "roundscale $w.upperfr.offsets.yoffset 2; popupcheckredraw visual_dialog_pop4"
        ttk::entry $w.upperfr.offsets.entx -width 4 -textvariable visoptions.xoffset -validate focus -takefocus 0 \
            -validatecommand "popupcheckredraw visual_dialog_pop3;my_validate %W 0 1 %P 2" \
            -invalidcommand "my_invalid %W;popupcheckredraw visual_dialog_pop3"
        ttk::entry $w.upperfr.offsets.enty -width 4 -textvariable visoptions.yoffset -validate focus -takefocus 0 \
            -validatecommand "popupcheckredraw visual_dialog_pop4;my_validate %W 0 1 %P 2" \
            -invalidcommand "my_invalid %W;popupcheckredraw visual_dialog_pop4"

	# pack $w.showclipsolution 
        
        pack $w.upperfr.size.sc $w.upperfr.size.ent -padx 4  -pady 12 -side left
        grid $w.upperfr.offsets.xoffset $w.upperfr.offsets.entx $w.upperfr.offsets.xlab -sticky nw -padx 4
        grid $w.upperfr.offsets.yoffset $w.upperfr.offsets.enty $w.upperfr.offsets.ylab -sticky nw -padx 4
	grid $w.upperfr.size $w.upperfr.offsets -sticky nw -pady 7 -padx 10
        grid anchor $w.upperfr center
	



	#	pack $w.lineartexture $w.numcols 



	ttk::labelframe $w.deform -relief groove -borderwidth 3 -text "Deformation settings"
	ttk::checkbutton $w.deform.cb -text "Deformation" \
	    -variable visoptions.deformation \
	    -command { Ng_Vis_Set parameters; redraw }

        ttk::label $w.deform.l -text "Scale: "
        ttk::spinbox $w.deform.sc1 -from 0 -to 1e99 -textvariable visoptions.scaledeform1 -width 5 \
            -command { Ng_Vis_Set parameters; redraw } \
            -validate focusout -validatecommand { Ng_Vis_Set parameters; redraw; string is double %P } \
            -invalidcommand { puts "invalid value, %P %s"; set visoptions.scaledeform1 1; }
        
	# tixControl $w.deform.sc1 -label "Scale: " -integer false \
	#     -variable visoptions.scaledeform1 \
	#     -command { Ng_Vis_Set parameters; redraw } \
	#     -options {
	# 	entry.width 6
	# 	label.width 7
	# 	label.anchor e
	#     }	

	ttk::scale $w.deform.sc2 -orient horizontal -length 100 -from 0 -to 1 \
	    -variable  visoptions.scaledeform2 \
	    -command { popupcheckredraw visual_dialog_pop5 }

	pack $w.deform -fill x -ipady 2 -pady 4 -ipady 3
	grid $w.deform.cb $w.deform.l $w.deform.sc1 $w.deform.sc2 -sticky nw -padx 4;# -side left -expand yes -anchor center -padx 4
	grid anchor $w.deform center
        grid columnconfigure $w.deform 0 -pad 20
        grid columnconfigure $w.deform 2 -pad 20

	ttk::labelframe $w.as -relief groove -borderwidth 3 -text "Scaling options"     
	ttk::checkbutton $w.as.autoscale -text "Autoscale" \
	    -variable visoptions.autoscale \
	    -command { Ng_Vis_Set parameters; redraw }

        ttk::label $w.as.lmin -text "Min-value"
        ttk::spinbox $w.as.smin -textvariable visoptions.mminval -width 5 -validate focus \
            -validatecommand "my_validatespinbox %W %P 10" \
            -command "Ng_Vis_Set parameters; redraw;" \
            -invalidcommand "my_invalidspinbox %W" -from -1e10 -to 1e10 -increment 0.001
            
        ttk::label $w.as.lmax -text "Max-value"
        ttk::spinbox $w.as.smax -textvariable visoptions.mmaxval -width 5 -validate focus \
            -validatecommand "Ng_Vis_Set parameters; redraw;my_validatespinbox %W %P 10" \
            -command "Ng_Vis_Set parameters; redraw;" \
            -invalidcommand "my_invalidspinbox %W;Ng_Vis_Set parameters; redraw" -from -1e10 -to 1e10 -increment 0.001
            
        #tixControl $w.as.minval -label "Min-value: " -integer false \
	    -variable visoptions.mminval \
	    -command { Ng_Vis_Set parametersrange; redraw } \
	    -options {
	#	entry.width 6
	#	label.width 12
	#	label.anchor e
	#    }	
	#tixControl $w.as.maxval -label "Max-value: " -integer false \
	    -variable visoptions.mmaxval \
	    -command { Ng_Vis_Set parametersrange; redraw } \
	    -options {
	#	entry.width 6
	#	label.width 12
	#	label.anchor e
	#    }	

	pack $w.as -fill x -pady 5 -ipady 3
	grid $w.as.autoscale $w.as.lmin $w.as.smin $w.as.lmax $w.as.smax -sticky nw -padx 4
        grid columnconfigure $w.as 0 -pad 20
        grid columnconfigure $w.as 2 -pad 20
        grid anchor $w.as center 




	ttk::frame $w.iso; #-relief groove -borderwidth 0
	pack $w.iso -anchor center;# -ipady 3

	ttk::labelframe $w.iso.cb -relief groove -borderwidth 3 -text "Iso lines / surfaces"
	pack $w.iso.cb -side left -pady 7 -fill y

	ttk::checkbutton $w.iso.cb.isolines -text "Iso-lines" \
	    -variable visoptions.isolines \
	    -command { Ng_Vis_Set parameters; redraw }
	#pack $w.iso.cb.isolines -side top -anchor w

	ttk::checkbutton $w.iso.cb.isosurf -text "Iso-Surface" \
	    -variable visoptions.isosurf \
	    -command { Ng_Vis_Set parameters; redraw }
	#pack $w.iso.cb.isosurf -side top -anchor w

        ttk::label $w.iso.cb.numisol -text "amount"
	ttk::scale $w.iso.cb.numiso -orient horizontal -length 100 -from 2 -to 50 \
	    -variable  visoptions.numiso \
	    -command "roundscale $w.iso.cb.numiso 0;popupcheckredraw visual_dialog_pop6"
        ttk::entry $w.iso.cb.entry -textvariable visoptions.numiso -width 3 \
            -validate focus -validatecommand "popupcheckredraw visual_dialog_pop6;\
            my_validate %W [$w.iso.cb.numiso cget -from] [$w.iso.cb.numiso cget -to] %P 0" \
            -invalidcommand "my_invalid %W;popupcheckredraw visual_dialog_pop6"
	    # -resolution 1    \
	    # -label "" \
        
        grid $w.iso.cb.isolines $w.iso.cb.numisol $w.iso.cb.entry -sticky nw -padx 4
        grid $w.iso.cb.isosurf  -sticky nw -padx 4
        grid $w.iso.cb.numiso  -sticky nw -padx 4 -columnspan 2 -column 1 -row 1 
	#pack $w.iso.cb.numisol $w.iso.cb.numiso  -anchor n


# 	scale $w.iso.subdiv -orient horizontal -length 100 -from 0 -to 5 \
# 	    -label "subdivision" \
# 	    -resolution 1    \
# 	    -variable  visoptions.subdivisions \
# 	    -command { popupcheckredraw visual_dialog_pop7  }
# #	    -command { puts "subdiv-vis"; Ng_Vis_Set parameters; puts "cal redraw"; redraw  }

	ttk::labelframe $w.iso.subdiv -text "Subdivision" -relief groove -borderwidth 3
	ttk::radiobutton $w.iso.subdiv.zero -text "0" -variable visoptions.subdivisions -value 0 \
	    -command { 
		#set visoptions.subdivisions  1; 
		Ng_Vis_Set parameters; redraw;
	}
	ttk::radiobutton $w.iso.subdiv.one -text "1" -variable visoptions.subdivisions -value 1 \
	    -command { 
		#set visoptions.subdivisions  1; 
		Ng_Vis_Set parameters; redraw;
	}
	ttk::radiobutton $w.iso.subdiv.two -text "2" -variable visoptions.subdivisions -value 2 \
	    -command { 
		#set visoptions.subdivisions  2; 
		Ng_Vis_Set parameters; redraw;
	}
	ttk::radiobutton $w.iso.subdiv.three -text "3" -variable visoptions.subdivisions -value 3 \
	    -command { 
		#set visoptions.subdivisions  3; 
		Ng_Vis_Set parameters; redraw;
	}
	ttk::radiobutton $w.iso.subdiv.four -text "4" -variable visoptions.subdivisions -value 4 \
	    -command { 
		#set visoptions.subdivisions  4; 
		Ng_Vis_Set parameters; redraw;
	}
	ttk::radiobutton $w.iso.subdiv.five -text "5" -variable visoptions.subdivisions -value 5 \
	    -command { 
		#set visoptions.subdivisions  5; 
		Ng_Vis_Set parameters; redraw;
	    }

 	ttk::label $w.iso.subdiv.text  -text "subdivision"

	pack $w.iso.subdiv -side right -fill y -padx 4 -pady 7 
# ; Ng_SetNextTimeStamp
	#pack $w.iso.subdiv.text  -side top
	#pack $w.iso.subdiv.zero  -side left
	#pack $w.iso.subdiv.one   -side left
	#pack $w.iso.subdiv.two   -side left
	#pack $w.iso.subdiv.three -side left
	#pack $w.iso.subdiv.four  -side left
	#pack $w.iso.subdiv.five  -side left
        grid $w.iso.subdiv.zero $w.iso.subdiv.one $w.iso.subdiv.two $w.iso.subdiv.three $w.iso.subdiv.four $w.iso.subdiv.five 
        grid anchor $w.iso.subdiv center
#	scale $w.iso.zpos -orient horizontal -length 100 -from 0 -to 1 \
#	    -label "z-position" \
#	    -resolution 0.01 \
#	    -variable visoptions.zposition \
#	    -command {
#		catch {NGS_Set zpos ${visoptions.zposition};}
#		redraw }
#	pack $w.iso.zpos -side right



	ttk::labelframe $w.redraw -relief groove -borderwidth 3 -text "Auto-redraw"
	ttk::checkbutton $w.redraw.auto -text "Auto-redraw after (sec)" \
	    -variable visoptions.autoredraw 

	# tixControl $w.redraw.val -integer false \
	#     -variable visoptions.autoredrawtime \
	#     -options {
	# 	entry.width 6
	# 	label.width 0
	# 	label.anchor w
	#     }	
        ttk::spinbox $w.redraw.val -textvariable visoptions.autoredrawtime -from 0 -to 100 -width 3 
	pack $w.redraw -fill x -ipady 3 -pady 7
	grid $w.redraw.auto  $w.redraw.val -sticky nw
        grid anchor $w.redraw center


        ttk::labelframe $w.lowerframe -text "Additional viewing options" -relief groove -borderwidth 3
        pack $w.lowerframe -fill x 
        set w $w.lowerframe 

        #pack [frame $w.f] -fill x
        #pack [ttk::frame $w.f1] -expand yes
        ttk::frame $w.f1
        
        set f [ttk::frame $w.f1.clipsol] 
        pack $f -anchor e
        menu $f.m
        ttk::menubutton $f.b -menu $f.m -width 12
        ttk::label $f.l -text "Clipping Plane Sol: "
        
        global visoptions.clipsolution        
        set clipsollabs(none) "None"
        set clipsollabs(scal) "Scalar Function"
        set clipsollabs(vec) "Vector Function"
        foreach i { none scal vec } {
            set textval $clipsollabs($i)
            $f.m add command -label "$textval" -command \
                "$f.b configure -text \"$textval\" ; set visoptions.clipsolution $i ; Ng_Vis_Set parameters ; redraw ; puts \"call redraw\" "                
        }

        pack $f.b $f.l -side right
        $f.m invoke $clipsollabs(${visoptions.clipsolution})

        
        # pack [ttk::frame $w.f1.scalfun] -anchor e
        set f [ttk::frame $w.f1.scalfun] 
        pack $f -anchor e
        menu $f.m
        ttk::menubutton $f.b -menu $f.m -width 12
        ttk::label $f.l -text "Scalar Function: "
        set scalentries [list none None]

	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
            set fname [Ng_Vis_Field getfieldname $i]
            set fcomp [Ng_Vis_Field getfieldcomponents $i]
            if { $fcomp == 1 } {
                lappend scalentries $fname:1 $fname
            } {
                for { set j 1 } { $j <= $fcomp } { incr j } {
                    lappend scalentries $fname:$j "$fname ($j)"
                }
                lappend scalentries $fname:0 "func ($fname)"
            }
        }
        global visoptions.scalfunction
        foreach { name textval } $scalentries {
            $f.m add command -label "$textval" -command \
                "$f.b configure -text \"$textval\" ; set visoptions.scalfunction $name ; Ng_Vis_Set parameters ; redraw ; "
        }
        
        pack $f.b $f.l -side right
        foreach { name textval } $scalentries {
            if { ${visoptions.scalfunction} == $name } {
                $f.m invoke $textval
            }
        }
        



        set f [ttk::frame $w.f1.vecfun]
        pack $f -anchor e
        menu $f.m
        ttk::menubutton $f.b -menu $f.m -width 12
        ttk::label $f.l -text "Vector Function: "
        set vecentries [list none None]

	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    set iscomplex [Ng_Vis_Field iscomplex $i]
	    set sdim [Ng_Vis_Field getdimension]
	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    if { ($fcomp == $sdim) || ($fcomp == 3) } {
                lappend vecentries $fname $fname
	    } 
        }
        global visoptions.vecfunction
        foreach { name textval } $vecentries {
            $f.m add command -label "$textval" -command \
                "$f.b configure -text \"$textval\" ; set visoptions.vecfunction $name ; Ng_Vis_Set parameters ; redraw ; "
        }
        
        pack $f.b $f.l -side right
        foreach { name textval } $vecentries {
            if { ${visoptions.vecfunction} == $name } {
                $f.m invoke $textval
            }
        }
        


        set f [ttk::frame $w.f1.evaluate]
        pack $f -anchor e

        menu $f.m
        ttk::menubutton $f.b -menu $f.m -width 12
        ttk::label $f.l -text "Evaluate: "
        
        global visoptions.evaluate        
        set evallabs(abs) "| |"
        set evallabs(abstens) "|tensor|"
        set evallabs(mises) "Mises"
        set evallabs(main) "Main"
        foreach i { abs abstens mises main } {
            set textval $evallabs($i)
            $f.m add command -label "$textval" -command \
                "$f.b configure -text \"$textval\" ; set visoptions.evaluate $i ; "
        }
        pack $f.b $f.l -side right
        $f.m invoke $evallabs(${visoptions.evaluate})

        
        
                
        pack [ttk::frame $w.f1.multidim] -fill x
        set f [ttk::frame $w.f1.multidim.f]
        pack $f -anchor e
        ttk::label $f.l1 -text "multidim-component: "
        ttk::spinbox $f.sb1 -from 0 -to 1e99 -textvariable visoptions.multidimcomponent -width 3 \
            -command { Ng_Vis_Set parameters; redraw }
        pack $f.l1 $f.sb1 -side left

        ttk::frame $w.fcb
        # the 2 main frames
        grid $w.f1 $w.fcb -sticky nw -padx 7 -ipady 3
        grid anchor $w center
        
        #pack $w.fcb 
        ttk::frame $w.fcb.cb
        pack $w.fcb.cb
        
	ttk::checkbutton $w.fcb.cb.showsurfsolution -text "Draw Surface Vectors" \
	    -variable visoptions.showsurfacesolution \
	    -command { Ng_Vis_Set parameters; redraw }

	ttk::checkbutton $w.fcb.cb.showcurves -text "Show Curves" \
	    -variable visoptions.drawpointcurves \
	    -command { Ng_Vis_Set parameters; redraw }

	ttk::checkbutton $w.fcb.cb.imaginary -text "Imaginary Part" \
	    -variable visoptions.imaginary \
	    -command { Ng_Vis_Set parameters; redraw }

	ttk::checkbutton $w.fcb.cb.logscale -text "Log Scale" \
	    -variable visoptions.logscale \
	    -command { Ng_Vis_Set parameters; redraw }

	ttk::checkbutton $w.fcb.cb.invcolor -text "Inverse Color" \
	    -variable visoptions.invcolor \
	    -command { Ng_Vis_Set parametersrange; redraw }


	ttk::frame $w.fcb.cb.texframe

	ttk::checkbutton $w.fcb.cb.texframe.usetexture -text "Use Textures (" \
	    -variable visoptions.usetexture \
	    -command { Ng_Vis_Set parameters; redraw }
	
	ttk::checkbutton $w.fcb.cb.texframe.lintexture -text "Linear )" \
	    -variable visoptions.lineartexture \
	    -command { Ng_Vis_Set parametersrange; redraw }
	
	ttk::checkbutton $w.fcb.cb.lineartexture -text "Use Linear Texture" \
	    -variable visoptions.lineartexture \
	    -command { Ng_Vis_Set parameters; redraw }
	
	scale $w.numcols -orient horizontal -length 100 -from 0 -to 50 \
	    -resolution 1   \
	    -variable  visoptions.numtexturecols \
	    -command { popupcheckredraw visual_dialog_pop1 }

	ttk::checkbutton $w.fcb.cb.showclipsolution -text "Draw Clipping Plane Solution" \
	    -variable visoptions.showclipsolution \
	    -command { Ng_Vis_Set parameters; redraw }


	ttk::checkbutton $w.fcb.cb.redrawperiodic -text "Animate periodic" \
	    -variable visoptions.redrawperiodic \
	    -command { 
		redrawperiodic
		Ng_Vis_Set parameters; 
		redraw 
	    }

	#pack $w.fcb.cb.showsurfsolution $w.fcb.cb.showcurves -anchor w
	#pack $w.fcb.cb.imaginary $w.fcb.cb.logscale $w.fcb.cb.texframe $w.fcb.cb.invcolor $w.fcb.cb.redrawperiodic -side top -anchor w
	#pack $w.fcb.cb.texframe.usetexture $w.fcb.cb.texframe.lintexture -side left -expand yes
	grid $w.fcb.cb.showsurfsolution -sticky nw
        grid $w.fcb.cb.showcurves -sticky nw
	grid $w.fcb.cb.imaginary -sticky nw
        grid $w.fcb.cb.logscale -sticky nw 
        grid $w.fcb.cb.texframe -sticky nw
        grid $w.fcb.cb.invcolor -sticky nw
        grid $w.fcb.cb.redrawperiodic -sticky nw
	pack $w.fcb.cb.texframe.usetexture $w.fcb.cb.texframe.lintexture -side left -expand yes
	


        set w .visoptions_dlg.main
	ttk::frame $w.bu;# -relief groove -borderwidth 3
	pack $w.bu  -pady 5 -padx 4
        
	ttk::button $w.bu.showsol -text "Show Solution" -command { 
	    set selectvisual solution
	    Ng_SetVisParameters
	    redraw
	}
	ttk::button $w.bu.clipping -text "Clipping" -command { 
	    clippingdialog; 
	}
	ttk::button $w.bu.fieldlines -text "Fieldlines" -command { 
	    fieldlinesdialog; 
	}

	ttk::button $w.bu.lineplot -text "2D Lineplot" -command {
	    lineplotdialog;
	}

	ttk::button $w.bu.done -text "Close" -command { 
	    destroy .visoptions_dlg
	}

	pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes 
        set w .visoptions_dlg

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Visualization"

    }
}



# proc reset_visual_dialog { } {
    
    # set w .visoptions_dlg
    
      # if {[winfo exists .visoptions_dlg] == 1} {
    
    
	  # destroy $w.scalfun $w.vecfun $w.evaluate $w.multidimcomp
	  # destroy $w.imaginary $w.logscale $w.texframe.usetexture $w.texframe.lintexture
          # destroy $w.texframe
          # destroy $w.invcolor $w.redrawperiodic
	  # destroy $w.bu  -pady 5
	  # destroy $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes
	  
	  
	  # checkbutton $w.imaginary -text "Imaginary Part" \
	      # -variable visoptions.imaginary \
	      # -command { Ng_Vis_Set parameters; redraw }

	  # frame $w.texframe

	  # checkbutton $w.texframe.usetexture -text "Use Textures (" \
	      # -variable visoptions.usetexture \
	      # -command { Ng_Vis_Set parameters; redraw }

	  # checkbutton $w.texframe.lintexture -text "Linear )" \
	      # -variable visoptions.lineartexture \
	      # -command { Ng_Vis_Set parameters; redraw }



	  
	  
	  # checkbutton $w.invcolor -text "Inverse Color" \
	      # -variable visoptions.invcolor \
	      # -command { Ng_Vis_Set parameters; redraw }

	  # checkbutton $w.logscale -text "Log Scale" \
	      # -variable visoptions.logscale \
	      # -command { Ng_Vis_Set parameters; redraw }
	  

	  
	  # checkbutton $w.redrawperiodic -text "Animate periodic" \
	      # -variable visoptions.redrawperiodic \
	      # -command { 
		  # redrawperiodic
		  # Ng_Vis_Set parameters; 
		  # redraw 
	      # }
	  

	# tixOptionMenu $w.scalfun -label "Scalar Function: " \
	    # -options {
		# label.width  18
		# label.anchor e
		# menubutton.width 12
	    # }

	# tixOptionMenu $w.vecfun -label "Vector Function: " \
	    # -options {
		# label.width  18
		# label.anchor e
		# menubutton.width 12
	    # }



	# $w.scalfun add command none -label None
	# for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    # set fname [Ng_Vis_Field getfieldname $i]
	    # set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    # if { $fcomp == 1 } {
		# $w.scalfun add command $fname.1 -label $fname
	    # } {
		# for { set j 1 } { $j <= $fcomp } { incr j } {
		    # $w.scalfun add command $fname.$j -label "$fname ($j)"
		# }
		# $w.scalfun add command $fname.0 -label "func ($fname)"
	    # }
	# }

	# $w.vecfun add command none -label None 
	# for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    # set fname [Ng_Vis_Field getfieldname $i]
	    # set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    # set iscomplex [Ng_Vis_Field iscomplex $i]
	    # set sdim [Ng_Vis_Field getdimension]
	    # if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    # if { ($fcomp == $sdim) || ($fcomp == 3) } {
		# $w.vecfun add command $fname -label $fname
	    # } 
	# }



	# $w.scalfun configure -variable visoptions.scalfunction 
	# $w.scalfun configure -command { Ng_Vis_Set parameters; redraw }
	# $w.vecfun configure -variable visoptions.vecfunction
	# $w.vecfun configure -command { Ng_Vis_Set parameters; redraw }


# #	puts "sclfunction = ${visoptions.scalfunction}"


	# tixOptionMenu $w.evaluate -label "Evaluate: " \
	    # -options {
		# label.width  18
		# label.anchor e
		# menubutton.width 12
	    # }	
	# $w.evaluate add command abs -label "|.|"
	# $w.evaluate add command abstens -label "|tensor|"
	# $w.evaluate add command mises -label "Mises"
	# $w.evaluate add command main  -label "Main"
	# $w.evaluate configure -variable visoptions.evaluate
	# $w.evaluate configure -command { 
	    # Ng_Vis_Set parameters; 
	    # redraw 
	# }

	# pack $w.scalfun $w.vecfun $w.evaluate

	# tixControl $w.multidimcomp -label "multidim-component: " -integer true \
	    # -variable visoptions.multidimcomponent -min 0 \
	    # -command { Ng_Vis_Set parameters; redraw } \
	    # -options {
		# entry.width 6
		# label.width 18
		# label.anchor e
	    # }	


          # pack $w.multidimcomp

          # pack $w.imaginary $w.logscale $w.texframe $w.invcolor $w.redrawperiodic
          # pack $w.texframe.usetexture $w.texframe.lintexture  -side left -expand yes


	# frame $w.bu
	# pack $w.bu  -pady 5

	# button $w.bu.showsol -text "Show Solution" -command { 
	    # set selectvisual solution
	    # Ng_SetVisParameters
	    # redraw
	# }
	# button $w.bu.clipping -text "Clipping" -command { 
	    # clippingdialog; 
	# }
	# button $w.bu.fieldlines -text "Fieldlines" -command { 
	    # fieldlinesdialog; 
	# }

	# button $w.bu.lineplot -text "2D Lineplot" -command {
	    # lineplotdialog;
	# }

	# button $w.bu.done -text "Close" -command { 
	    # destroy .visoptions_dlg
	# }

	# pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes

	# wm withdraw $w
	# wm deiconify $w


      # }

# }



