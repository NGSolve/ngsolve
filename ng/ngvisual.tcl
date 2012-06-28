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
	
	frame $w.filesettings -relief  groove -borderwidth 3
	frame $w.filesettings.title
	radiobutton $w.filesettings.title.choose -variable visoptions.lineplotsource \
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
	

	tixOptionMenu $w.filesettings.latestevals -label "Use Evaluate Results: " \
	    -options {
		label.width  25
		label.anchor e
		menubutton.width 40
	    } 
	
	for {set i 0} {$i < [llength ${visoptions.evaluatefilenames}]} {incr i} {
	    $w.filesettings.latestevals add command $i \
		-label "[lindex ${visoptions.evaluatefiledescriptions} $i] ([lindex ${visoptions.evaluatefilenames} $i])"
	}
	$w.filesettings.latestevals config -variable visoptions.lineplotselectedeval

	pack $w.filesettings.latestevals
	
	frame $w.filesettings.sfn

	button $w.filesettings.sfn.bb -text "Browse" \
	    -command { set visoptions.lineplotfile [tk_getOpenFile] }

	
	entry $w.filesettings.sfn.fn -width 50 -relief sunken \
	    -textvariable visoptions.lineplotfile

	pack $w.filesettings.sfn.bb $w.filesettings.sfn.fn -side left

	pack $w.filesettings.sfn

	button $w.filesettings.refresh -text "Refresh" -command {
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


	frame $w.filesettings.using

	global visoptions.lineplotdatadescr
	
	tixOptionMenu $w.filesettings.using.xco -label "X-Coord:"\
	    -options {
		label.width  8
		label.anchor e
		menubutton.width 15
	    } 
	for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
	    $w.filesettings.using.xco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	}
	$w.filesettings.using.xco config -variable visoptions.lineplotusingx
	
	tixOptionMenu $w.filesettings.using.yco -label "Y-Coord:"\
	    -options {
		label.width  8
		label.anchor e
		menubutton.width 15
	    } 
	for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {
	    $w.filesettings.using.yco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]
	}
	$w.filesettings.using.yco config -variable visoptions.lineplotusingy

	global visoptions.lineplotxcoordselector
	global visoptions.lineplotycoordselector
	set visoptions.lineplotxcoordselector $w.filesettings.using.xco
	set visoptions.lineplotycoordselector $w.filesettings.using.yco
	

	pack $w.filesettings.using.xco $w.filesettings.using.yco -side left
	pack $w.filesettings.using

	pack $w.filesettings -fill x -ipady 3
	
	frame $w.settings -relief  groove -borderwidth 3
	label $w.settings.title -text "\nSettings\n"
	pack $w.settings.title

	frame $w.settings.minmax 
	checkbutton $w.settings.minmax.autoscale -text "Autoscale" -variable visoptions.lineplotautoscale
	tixControl $w.settings.minmax.xmin -label "Min. x: " \
	    -integer false -variable visoptions.lineplotxmin \
	    -options {
		entry.width 6
		label.width 8
		label.anchor e
	    }	
	tixControl $w.settings.minmax.xmax -label "Max. x: " \
	    -integer false -variable visoptions.lineplotxmax \
	    -options {
		entry.width 6
		label.width 8
		label.anchor e
	    }	
	tixControl $w.settings.minmax.ymin -label "Min. y: " \
	    -integer false -variable visoptions.lineplotymin \
	    -options {
		entry.width 6
		label.width 8
		label.anchor e
	    }	
	tixControl $w.settings.minmax.ymax -label "Max. y: " \
	    -integer false -variable visoptions.lineplotymax \
	    -options {
		entry.width 6
		label.width 8
		label.anchor e
	    }	
	    

	pack $w.settings.minmax.autoscale $w.settings.minmax.xmin $w.settings.minmax.xmax \
	    $w.settings.minmax.ymin $w.settings.minmax.ymax -side left

	pack $w.settings.minmax

	
	label $w.settings.empty1 -text ""
	pack $w.settings.empty1

	frame $w.settings.plotsize

	tixControl $w.settings.plotsize.xsize -label "Plotsize  x: "\
	    -integer true -variable visoptions.lineplotsizex \
	    -options {
		entry.width 6
		label.width 13
		label.anchor e
	    }	
	tixControl $w.settings.plotsize.ysize -label "y: "\
	    -integer true -variable visoptions.lineplotsizey \
	    -options {
		entry.width 6
		label.width 3
		label.anchor e
	    }	

	pack $w.settings.plotsize.xsize $w.settings.plotsize.ysize -side left

	pack $w.settings.plotsize

	label $w.settings.empty2 -text ""
	pack $w.settings.empty2
	

	tixOptionMenu $w.settings.color -label "Linecolor: " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 
	foreach step { red black blue green yellow } {
	    $w.settings.color add command $step -label $step
	}
	$w.settings.color config -variable visoptions.lineplotcolor

	pack $w.settings.color


	pack $w.settings

	set datax ""
	set datay ""
	set xmin 0
	set xmax 0
	set ymin 0
	set ymax 0

	frame $w.plots -relief  groove -borderwidth 3

	tixOptionMenu $w.plots.selplot -label "Selected Plot: " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 
	$w.plots.selplot add command none -label "None"

	$w.plots.selplot config -variable visoptions.lineplotselected

	global visoptions.lineplotselector
	set visoptions.lineplotselector $w.plots.selplot

	

	button $w.plots.new -text "Generate New Plot" -command {
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
	
	button $w.plots.addto -text "Add to Selected Plot" -command {
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
	


	button $w.close -text "Close" -command "destroy $w"
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

	tixNoteBook $w.nb -ipadx 6 -ipady 6

	$w.nb add draw -label "Draw"
	$w.nb add settings -label "Settings"

	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top


	# Main Window

	set f [$w.nb subwidget draw]
	

	frame $f.general

	
	checkbutton $f.general.enable -text "Enable Fieldlines" \
	    -variable visoptions.drawfieldlines \
	    -command { 
		# set visoptions.redrawperiodic ${visoptions.drawfieldlines}
		# redrawperiodic
		# redrawperiodic # sonst 
		Ng_Vis_Set parameters; 
		redraw 
	    }
	tixControl $f.general.num -label "Num: " -integer true \
	    -variable visoptions.numfieldlines \
	    -command { Ng_Vis_Set parameters; redraw } \
	    -options {
		entry.width 6
		label.width 12
		label.anchor e
	    }	

	
	pack $f.general.enable $f.general.num -side left

	pack $f.general

	label $f.labe0 -text " "

	pack $f.labe0

	frame $f.general1
	
	checkbutton $f.general1.randomstart -text "Field dependent density    " \
	    -variable visoptions.fieldlinesrandomstart \
	    -command { Ng_Vis_Set parameters; redraw}

	
	checkbutton $f.general1.redrawperiodic -text "Animate periodic" \
	    -variable visoptions.redrawperiodic \
	    -command { 
		redrawperiodic
		Ng_Vis_Set parameters; 
		redraw 
	    }

	pack $f.general1.randomstart $f.general1.redrawperiodic -side left

	pack $f.general1



	label $f.lab0 -text " "

	pack $f.lab0


	
	tixOptionMenu $f.vecfun -label "Vector Function: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }
	$f.vecfun add command none -label None 
	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    set iscomplex [Ng_Vis_Field iscomplex $i]
	    set sdim [Ng_Vis_Field getdimension]
	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    if { ($fcomp == $sdim) || ($fcomp == 3) } {
		$f.vecfun add command $fname -label $fname
	    } 
	}
	$f.vecfun configure -variable visoptions.fieldlinesvecfunction
	$f.vecfun configure -command { Ng_Vis_Set parameters; redraw }

	pack $f.vecfun
	


	label $f.lab00 -text " "

	pack $f.lab00

	frame $f.phasesettings

	checkbutton $f.phasesettings.onephase -text "Fix Phase" -variable visoptions.fieldlinesonlyonephase
	scale $f.phasesettings.phase -orient horizontal -length 300 -from 0 -to 360 \
	    -label "phi" \
	    -resolution 1 \
	    -variable visoptions.fieldlinesphase \
	    -command { popupcheckredraw3 fieldlinesdialog_pop1 }

	pack $f.phasesettings.onephase $f.phasesettings.phase -side left
	
	pack $f.phasesettings



	label $f.lab1 -text " "

	pack $f.lab1


	
	frame $f.boxsettings -relief groove -borderwidth 3
	frame $f.boxsettings.title
	radiobutton $f.boxsettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value box -text "Startpoints in Box"

	pack $f.boxsettings.title.choose -side left

	pack $f.boxsettings.title

	frame $f.boxsettings.points

	label $f.boxsettings.points.lab2 -text "Pmin";
	entry $f.boxsettings.points.ent1x -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap1x
	entry $f.boxsettings.points.ent1y -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap1y
	entry $f.boxsettings.points.ent1z -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap1z
	label $f.boxsettings.points.lab3 -text "   Pmax";
	entry $f.boxsettings.points.ent2x -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap2x
	entry $f.boxsettings.points.ent2y -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap2y
	entry $f.boxsettings.points.ent2z -width 8 -relief sunken \
	    -textvariable visoptions.fieldlinesstartareap2z
	
	pack $f.boxsettings.points
	pack $f.boxsettings.points.lab2 $f.boxsettings.points.ent1x $f.boxsettings.points.ent1y $f.boxsettings.points.ent1z -side left
	pack $f.boxsettings.points.lab3 $f.boxsettings.points.ent2x $f.boxsettings.points.ent2y $f.boxsettings.points.ent2z -side left

	button $f.boxsettings.settobb -text "Bounding Box" -command {
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


	frame $f.facesettings -relief groove -borderwidth 3
	frame $f.facesettings.title
	radiobutton $f.facesettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value face -text "Startpoints on Face"

	pack $f.facesettings.title.choose -side left

	pack $f.facesettings.title
	
	frame $f.facesettings.index
	label $f.facesettings.index.lab -text "face index:"
	label $f.facesettings.index.ent -text 1 -padx 4

	pack $f.facesettings.index.lab $f.facesettings.index.ent -side left

	pack $f.facesettings.index
	
	pack $f.facesettings -fill x -ipady 3


	global visoptions.fieldlinesfilename

	frame $f.filesettings -relief  groove -borderwidth 3
	frame $f.filesettings.title
	radiobutton $f.filesettings.title.choose -variable visoptions.fieldlinesstartarea \
	    -value file -text "Startpoints from File"

	pack $f.filesettings.title.choose -side left

	pack $f.filesettings.title

	frame $f.filesettings.sfn

	button $f.filesettings.sfn.bb -text "Browse" \
	    -command {
		set types {
		    { "Netgen Fieldlines" {.nef} }
		}
		set visoptions.fieldlinesfilename [tk_getOpenFile -filetypes $types -defaultextension ".nef"]
	    }

	
	entry $f.filesettings.sfn.fn -width 50 -relief sunken \
	    -textvariable visoptions.fieldlinesfilename

	pack $f.filesettings.sfn.bb $f.filesettings.sfn.fn -side left

	pack $f.filesettings.sfn

	pack $f.filesettings -fill x -ipady 3
	


	
	# Settings

	set g [$w.nb subwidget settings]

	frame $g.linesettings -relief groove -borderwidth 3
	label $g.linesettings.title -text "\nLine Settings\n"
	tixControl $g.linesettings.length -label "rel. Length: " -integer false \
	    -variable visoptions.fieldlineslength -min 0.00001 -max 10000 -step 0.1 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }

	tixControl $g.linesettings.maxpoints -label "max. Points: " -integer true \
	    -variable visoptions.fieldlinesmaxpoints -min 0 -max 10000 -step 1 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }

	tixControl $g.linesettings.thick -label "rel. Thickness: " -integer false \
	    -variable visoptions.fieldlinesthickness -min 1e-10 -max 0.5 -step 0.001 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }

	pack $g.linesettings.title $g.linesettings.length $g.linesettings.maxpoints $g.linesettings.thick

	pack $g.linesettings -fill x -ipady 3


	


	global visoptions.fieldlinestolerance

	frame $g.odesettings -relief groove -borderwidth 3
	label $g.odesettings.title -text "\nODE Settings\n"
	tixControl $g.odesettings.tol -label "rel. Tolerance: " -integer false \
	    -variable visoptions.fieldlinestolerance -min 0.00001 -max 1 -step 0.01 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	


	tixOptionMenu $g.odesettings.rktype -label "RK-Type " \
	    -options {
		label.width  20
		label.anchor e
		menubutton.width 25
	    }
	$g.odesettings.rktype add command euler -label "Euler, order 1"
	$g.odesettings.rktype add command eulercauchy -label "Euler-Cauchy, order 2"
	$g.odesettings.rktype add command simpson -label "Simpson, order 3"
	$g.odesettings.rktype add command crungekutta -label "classical Runge-Kutta, order 4"
	$g.odesettings.rktype configure -variable visoptions.fieldlinesrktype
	$g.odesettings.rktype configure -command { Ng_Vis_Set parameters; redraw }
	
	pack $g.odesettings.title $g.odesettings.tol $g.odesettings.rktype

	pack $g.odesettings -fill x -ipady 3



	# buttons
	

	frame $w.bu
	pack $w.bu 

	button $w.bu.calc -text "Build Fieldlines" -command { 
	    if { ${visoptions.fieldlinesvecfunction} == "none" } {
		bgerror "Please select the vector function first!"
	    } {
		set visoptions.drawfieldlines 1
		Ng_Vis_Set parameters
		Ng_BuildFieldLines
		redraw 
	    }
	}

	button $w.bu.help -text "Help" -command {
	    if {[winfo exists .fieldlines_help] == 1} {
		wm withdraw .fieldlines_help
		wm deiconify .fieldlines_help
		focus .fieldlines_help
	    } {
		toplevel .fieldlines_help

		tixScrolledText .fieldlines_help.ht -scrollbar y
		set text [.fieldlines_help.ht subwidget text]

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

	button $w.bu.cancel -text "Done" -command "destroy $w"
	pack $w.bu.calc $w.bu.help $w.bu.cancel -side left -expand yes
	
	
	wm withdraw $w
	wm geom $w +200+100
	wm deiconify $w
	wm title $w "Fieldlines"
	#    grab $w
	focus $w

    }

    global visoptions.fieldlinesstartface

    
    set f [$w.nb subwidget draw]
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


	


	frame $w.grid -relief groove -borderwidth 3
	# change to: max gridsize 200
	scale $w.grid.size -orient horizontal -length 100 -from 1 -to 200 \
	    -label "Grid" \
	    -resolution 1    \
	    -variable  visoptions.gridsize \
	    -command { popupcheckredraw visual_dialog_pop2 }
    

	# x- and y- offset
	scale $w.grid.xoffset -orient horizontal -length 80 -from 0 -to 1 \
	    -label "x-Offset" \
	    -resolution 0.05    \
	    -variable  visoptions.xoffset \
	    -command { popupcheckredraw visual_dialog_pop3 }

	scale $w.grid.yoffset -orient horizontal -length 80 -from 0 -to 1 \
	    -label "y-Offset" \
	    -resolution 0.05    \
	    -variable  visoptions.yoffset \
	    -command { popupcheckredraw visual_dialog_pop4 }


	# pack $w.showclipsolution 
	pack $w.grid -fill x -ipady 3
	pack $w.grid.size $w.grid.xoffset $w.grid.yoffset -side left -expand yes



	#	pack $w.lineartexture $w.numcols 



	frame $w.deform -relief groove -borderwidth 3
	checkbutton $w.deform.cb -text "Deformation" \
	    -variable visoptions.deformation \
	    -command { Ng_Vis_Set parameters; redraw }

	tixControl $w.deform.sc1 -label "Scale: " -integer false \
	    -variable visoptions.scaledeform1 \
	    -command { Ng_Vis_Set parameters; redraw } \
	    -options {
		entry.width 6
		label.width 7
		label.anchor e
	    }	

	scale $w.deform.sc2 -orient horizontal -length 100 -from 0 -to 1 \
	    -resolution 0.01    \
	    -variable  visoptions.scaledeform2 \
	    -command { popupcheckredraw visual_dialog_pop5 }

	pack $w.deform -fill x -ipady 2
	pack $w.deform.cb $w.deform.sc1 $w.deform.sc2 -side left -expand yes
	

	frame $w.as -relief groove -borderwidth 3
	checkbutton $w.as.autoscale -text "Autoscale" \
	    -variable visoptions.autoscale \
	    -command { Ng_Vis_Set parameters; redraw }

	tixControl $w.as.minval -label "Min-value: " -integer false \
	    -variable visoptions.mminval \
	    -command { Ng_Vis_Set parametersrange; redraw } \
	    -options {
		entry.width 6
		label.width 12
		label.anchor e
	    }	
	tixControl $w.as.maxval -label "Max-value: " -integer false \
	    -variable visoptions.mmaxval \
	    -command { Ng_Vis_Set parametersrange; redraw } \
	    -options {
		entry.width 6
		label.width 12
		label.anchor e
	    }	

	pack $w.as -fill x -ipady 3
	pack $w.as.autoscale $w.as.minval $w.as.maxval -side left




	frame $w.iso -relief groove -borderwidth 3
	pack $w.iso -fill x -ipady 3

	frame $w.iso.cb
	pack $w.iso.cb -side left

	checkbutton $w.iso.cb.isolines -text "Iso-lines" \
	    -variable visoptions.isolines \
	    -command { Ng_Vis_Set parameters; redraw }
	pack $w.iso.cb.isolines -side top

	checkbutton $w.iso.cb.isosurf -text "Iso-Surface" \
	    -variable visoptions.isosurf \
	    -command { Ng_Vis_Set parameters; redraw }
	pack $w.iso.cb.isosurf -side top



	scale $w.iso.numiso -orient horizontal -length 100 -from 2 -to 50 \
	    -label "" \
	    -resolution 1    \
	    -variable  visoptions.numiso \
	    -command { popupcheckredraw visual_dialog_pop6 }

	pack $w.iso.numiso -side left


# 	scale $w.iso.subdiv -orient horizontal -length 100 -from 0 -to 5 \
# 	    -label "subdivision" \
# 	    -resolution 1    \
# 	    -variable  visoptions.subdivisions \
# 	    -command { popupcheckredraw visual_dialog_pop7  }
# #	    -command { puts "subdiv-vis"; Ng_Vis_Set parameters; puts "cal redraw"; redraw  }

	frame $w.iso.subdiv
	radiobutton $w.iso.subdiv.zero -text "0" -variable visoptions.subdivisions -value 0 \
	    -command { 
		#set visoptions.subdivisions  1; 
		Ng_Vis_Set parameters; redraw;
	}
	radiobutton $w.iso.subdiv.one -text "1" -variable visoptions.subdivisions -value 1 \
	    -command { 
		#set visoptions.subdivisions  1; 
		Ng_Vis_Set parameters; redraw;
	}
	radiobutton $w.iso.subdiv.two -text "2" -variable visoptions.subdivisions -value 2 \
	    -command { 
		#set visoptions.subdivisions  2; 
		Ng_Vis_Set parameters; redraw;
	}
	radiobutton $w.iso.subdiv.three -text "3" -variable visoptions.subdivisions -value 3 \
	    -command { 
		#set visoptions.subdivisions  3; 
		Ng_Vis_Set parameters; redraw;
	}
	radiobutton $w.iso.subdiv.four -text "4" -variable visoptions.subdivisions -value 4 \
	    -command { 
		#set visoptions.subdivisions  4; 
		Ng_Vis_Set parameters; redraw;
	}
	radiobutton $w.iso.subdiv.five -text "5" -variable visoptions.subdivisions -value 5 \
	    -command { 
		#set visoptions.subdivisions  5; 
		Ng_Vis_Set parameters; redraw;
	    }

 	label $w.iso.subdiv.text  -text "subdivision"

	pack $w.iso.subdiv -side right -ipadx 10

# ; Ng_SetNextTimeStamp
	pack $w.iso.subdiv.text -side top
	pack $w.iso.subdiv.zero $w.iso.numiso -side left
	pack $w.iso.subdiv.one $w.iso.numiso -side left
	pack $w.iso.subdiv.two $w.iso.numiso -side left
	pack $w.iso.subdiv.three $w.iso.numiso -side left
	pack $w.iso.subdiv.four $w.iso.numiso -side left
	pack $w.iso.subdiv.five $w.iso.numiso -side left


#	scale $w.iso.zpos -orient horizontal -length 100 -from 0 -to 1 \
#	    -label "z-position" \
#	    -resolution 0.01 \
#	    -variable visoptions.zposition \
#	    -command {
#		catch {NGS_Set zpos ${visoptions.zposition};}
#		redraw }
#	pack $w.iso.zpos -side right



	frame $w.redraw -relief groove -borderwidth 3
	checkbutton $w.redraw.auto -text "Auto-redraw" \
	    -variable visoptions.autoredraw 

	tixControl $w.redraw.val -label " after (sec) " -integer false \
	    -variable visoptions.autoredrawtime \
	    -options {
		entry.width 6
		label.width 0
		label.anchor w
	    }	

	pack $w.redraw -fill x -ipady 3
	pack $w.redraw.auto  $w.redraw.val -side left



	tixControl $w.redraw.simtime -label " Simulation Time (1e-6 s)" -integer false \
	    -variable visoptions.simulationtime \
	    -command { 
		Ng_Vis_Set time ${visoptions.simulationtime}; 
		catch {NGS_Set time ${visoptions.simulationtime};}
		redraw } \
	    -options {
		entry.width 6
		label.width 0
		label.anchor w
	    }	
	pack $w.redraw.simtime -side left
	



	tixOptionMenu $w.clipsol -label "Clipping Plane Sol: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }

	set none 1
	$w.clipsol add command none -label None
	$w.clipsol add command scal -label "Scalar Function"
	$w.clipsol add command vec -label "Vector Function"

	$w.clipsol configure -variable visoptions.clipsolution
	$w.clipsol configure -command { Ng_Vis_Set parameters; redraw }

	pack $w.clipsol




	tixOptionMenu $w.scalfun -label "Scalar Function: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }

	tixOptionMenu $w.vecfun -label "Vector Function: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }


	$w.scalfun add command none -label None
	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    if { $fcomp == 1 } {
		$w.scalfun add command $fname:1 -label $fname
	    } {
		for { set j 1 } { $j <= $fcomp } { incr j } {
		    $w.scalfun add command $fname:$j -label "$fname ($j)"
		}
		$w.scalfun add command $fname:0 -label "func ($fname)"
	    }
	}

	$w.vecfun add command none -label None 
	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    set iscomplex [Ng_Vis_Field iscomplex $i]
	    set sdim [Ng_Vis_Field getdimension]
	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    if { ($fcomp == $sdim) || ($fcomp == 3) } {
		$w.vecfun add command $fname -label $fname
	    } 
	}

	$w.scalfun configure -variable visoptions.scalfunction 
	$w.scalfun configure -command { Ng_Vis_Set parameters; redraw }
	$w.vecfun configure -variable visoptions.vecfunction
	$w.vecfun configure -command { Ng_Vis_Set parameters; redraw }


	tixOptionMenu $w.evaluate -label "Evaluate: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }	
	$w.evaluate add command abs -label "|.|"
	$w.evaluate add command abstens -label "|tensor|"
	$w.evaluate add command mises -label "Mises"
	$w.evaluate add command main  -label "Main"
	$w.evaluate configure -variable visoptions.evaluate
	$w.evaluate configure -command { 
	    Ng_Vis_Set parameters; 
	    redraw 
	}

	pack $w.scalfun $w.vecfun $w.evaluate

	tixControl $w.multidimcomp -label "multidim-component: " -integer true \
	    -variable visoptions.multidimcomponent -min 0 \
	    -command { Ng_Vis_Set parameters; redraw } \
	    -options {
		entry.width 6
		label.width 18
		label.anchor e
	    }	

	pack $w.multidimcomp


	checkbutton $w.showsurfsolution -text "Draw Surface Vectors" \
	    -variable visoptions.showsurfacesolution \
	    -command { Ng_Vis_Set parameters; redraw }

	checkbutton $w.showcurves -text "Show Curves" \
	    -variable visoptions.drawpointcurves \
	    -command { Ng_Vis_Set parameters; redraw }

	checkbutton $w.imaginary -text "Imaginary Part" \
	    -variable visoptions.imaginary \
	    -command { Ng_Vis_Set parameters; redraw }

	checkbutton $w.logscale -text "Log Scale" \
	    -variable visoptions.logscale \
	    -command { Ng_Vis_Set parameters; redraw }

	checkbutton $w.invcolor -text "Inverse Color" \
	    -variable visoptions.invcolor \
	    -command { Ng_Vis_Set parametersrange; redraw }


	frame $w.texframe

	checkbutton $w.texframe.usetexture -text "Use Textures (" \
	    -variable visoptions.usetexture \
	    -command { Ng_Vis_Set parameters; redraw }
	
	checkbutton $w.texframe.lintexture -text "Linear )" \
	    -variable visoptions.lineartexture \
	    -command { Ng_Vis_Set parametersrange; redraw }


	
	checkbutton $w.lineartexture -text "Use Linear Texture" \
	    -variable visoptions.lineartexture \
	    -command { Ng_Vis_Set parameters; redraw }
	
	scale $w.numcols -orient horizontal -length 100 -from 0 -to 50 \
	    -resolution 1   \
	    -variable  visoptions.numtexturecols \
	    -command { popupcheckredraw visual_dialog_pop1 }

	checkbutton $w.showclipsolution -text "Draw Clipping Plane Solution" \
	    -variable visoptions.showclipsolution \
	    -command { Ng_Vis_Set parameters; redraw }


	checkbutton $w.redrawperiodic -text "Animate periodic" \
	    -variable visoptions.redrawperiodic \
	    -command { 
		redrawperiodic
		Ng_Vis_Set parameters; 
		redraw 
	    }


	pack $w.showsurfsolution $w.showcurves
	pack $w.imaginary $w.logscale $w.texframe $w.invcolor $w.redrawperiodic
	pack $w.texframe.usetexture $w.texframe.lintexture -side left -expand yes
	



	frame $w.bu
	pack $w.bu  -pady 5

	button $w.bu.showsol -text "Show Solution" -command { 
	    set selectvisual solution
	    Ng_SetVisParameters
	    redraw
	}
	button $w.bu.clipping -text "Clipping" -command { 
	    clippingdialog; 
	}
	button $w.bu.fieldlines -text "Fieldlines" -command { 
	    fieldlinesdialog; 
	}

	button $w.bu.lineplot -text "2D Lineplot" -command {
	    lineplotdialog;
	}

	button $w.bu.done -text "Close" -command { 
	    destroy .visoptions_dlg
	}

	pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Visualization"

    }
}



proc reset_visual_dialog { } {
    
    set w .visoptions_dlg
    
      if {[winfo exists .visoptions_dlg] == 1} {
    
    
	  destroy $w.scalfun $w.vecfun $w.evaluate $w.multidimcomp
	  destroy $w.imaginary $w.logscale $w.texframe.usetexture $w.texframe.lintexture
          destroy $w.texframe
          destroy $w.invcolor $w.redrawperiodic
	  destroy $w.bu  -pady 5
	  destroy $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes
	  
	  
	  checkbutton $w.imaginary -text "Imaginary Part" \
	      -variable visoptions.imaginary \
	      -command { Ng_Vis_Set parameters; redraw }

	  frame $w.texframe

	  checkbutton $w.texframe.usetexture -text "Use Textures (" \
	      -variable visoptions.usetexture \
	      -command { Ng_Vis_Set parameters; redraw }

	  checkbutton $w.texframe.lintexture -text "Linear )" \
	      -variable visoptions.lineartexture \
	      -command { Ng_Vis_Set parameters; redraw }



	  
	  
	  checkbutton $w.invcolor -text "Inverse Color" \
	      -variable visoptions.invcolor \
	      -command { Ng_Vis_Set parameters; redraw }

	  checkbutton $w.logscale -text "Log Scale" \
	      -variable visoptions.logscale \
	      -command { Ng_Vis_Set parameters; redraw }
	  

	  
	  checkbutton $w.redrawperiodic -text "Animate periodic" \
	      -variable visoptions.redrawperiodic \
	      -command { 
		  redrawperiodic
		  Ng_Vis_Set parameters; 
		  redraw 
	      }
	  

	tixOptionMenu $w.scalfun -label "Scalar Function: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }

	tixOptionMenu $w.vecfun -label "Vector Function: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }



	$w.scalfun add command none -label None
	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    if { $fcomp == 1 } {
		$w.scalfun add command $fname.1 -label $fname
	    } {
		for { set j 1 } { $j <= $fcomp } { incr j } {
		    $w.scalfun add command $fname.$j -label "$fname ($j)"
		}
		$w.scalfun add command $fname.0 -label "func ($fname)"
	    }
	}

	$w.vecfun add command none -label None 
	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {
	    set fname [Ng_Vis_Field getfieldname $i]
	    set fcomp [Ng_Vis_Field getfieldcomponents $i]
	    set iscomplex [Ng_Vis_Field iscomplex $i]
	    set sdim [Ng_Vis_Field getdimension]
	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }
	    if { ($fcomp == $sdim) || ($fcomp == 3) } {
		$w.vecfun add command $fname -label $fname
	    } 
	}



	$w.scalfun configure -variable visoptions.scalfunction 
	$w.scalfun configure -command { Ng_Vis_Set parameters; redraw }
	$w.vecfun configure -variable visoptions.vecfunction
	$w.vecfun configure -command { Ng_Vis_Set parameters; redraw }


#	puts "sclfunction = ${visoptions.scalfunction}"


	tixOptionMenu $w.evaluate -label "Evaluate: " \
	    -options {
		label.width  18
		label.anchor e
		menubutton.width 12
	    }	
	$w.evaluate add command abs -label "|.|"
	$w.evaluate add command abstens -label "|tensor|"
	$w.evaluate add command mises -label "Mises"
	$w.evaluate add command main  -label "Main"
	$w.evaluate configure -variable visoptions.evaluate
	$w.evaluate configure -command { 
	    Ng_Vis_Set parameters; 
	    redraw 
	}

	pack $w.scalfun $w.vecfun $w.evaluate

	tixControl $w.multidimcomp -label "multidim-component: " -integer true \
	    -variable visoptions.multidimcomponent -min 0 \
	    -command { Ng_Vis_Set parameters; redraw } \
	    -options {
		entry.width 6
		label.width 18
		label.anchor e
	    }	


          pack $w.multidimcomp

          pack $w.imaginary $w.logscale $w.texframe $w.invcolor $w.redrawperiodic
          pack $w.texframe.usetexture $w.texframe.lintexture  -side left -expand yes


	frame $w.bu
	pack $w.bu  -pady 5

	button $w.bu.showsol -text "Show Solution" -command { 
	    set selectvisual solution
	    Ng_SetVisParameters
	    redraw
	}
	button $w.bu.clipping -text "Clipping" -command { 
	    clippingdialog; 
	}
	button $w.bu.fieldlines -text "Fieldlines" -command { 
	    fieldlinesdialog; 
	}

	button $w.bu.lineplot -text "2D Lineplot" -command {
	    lineplotdialog;
	}

	button $w.bu.done -text "Close" -command { 
	    destroy .visoptions_dlg
	}

	pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes

	wm withdraw $w
	wm deiconify $w


      }

}



