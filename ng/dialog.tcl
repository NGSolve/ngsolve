proc meshingoptionsdialog { } {

    set w .options_dlg
    
    if {[winfo exists .options_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {

	toplevel $w

#	global options.meshsize

	tixNoteBook $w.nb -ipadx 6 -ipady 6
	
	$w.nb add general -label "General" -underline 0
	$w.nb add meshsize -label "Mesh Size"   -underline 0
	$w.nb add chartopt -label "STL Charts" -underline 0
	$w.nb add optimizer -label "Optimizer"   -underline 0
	$w.nb add insider -label "Insider"   -underline 0
	$w.nb add debug -label "Debug"   -underline 0


	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	


	# General meshing options

	set f [$w.nb subwidget general]

	set finevals { 1 2 3 4 5 6 }
	set finelabs(1) "very coarse" 
	set finelabs(2) "coarse" 
	set finelabs(3) "moderate" 
	set finelabs(4) "fine" 
	set finelabs(5) "very fine" 
	set finelabs(6) "user defined"

	tixOptionMenu $f.fine -label "Mesh granularity : " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 


	foreach finev $finevals {
	    $f.fine add command $finev -label $finelabs($finev)
	}
	$f.fine config -variable meshoptions.fineness
	$f.fine config -command { setgranularity }
	global meshoptions.fineness
#	setgranularity ${meshoptions.fineness}
	pack $f.fine


	


	set mgsteps { ag me ms os mv ov }
	set mgsteplabel(ag) "Analyze Geometry"
	set mgsteplabel(me) "Mesh Edges"
	set mgsteplabel(ms) "Mesh Surface"
	set mgsteplabel(os) "Optimize Surface"
	set mgsteplabel(mv) "Mesh Volume"
	set mgsteplabel(ov) "Optimize Volume"

	
	tixOptionMenu $f.first -label "First Step : " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 

	tixOptionMenu $f.last -label "Last Step : " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 

	foreach step $mgsteps {
	    $f.first add command $step -label $mgsteplabel($step)
	    $f.last add command $step -label $mgsteplabel($step)
	}

	$f.first config -variable meshoptions.firststep 
	$f.last config  -variable meshoptions.laststep 

	pack $f.first $f.last
	




	set msg(0) "None"
	set msg(1) "Least"
	set msg(2) "Little"
	set msg(3) "Moderate"
	set msg(4) "Much"
	set msg(5) "Most"
	
	tixOptionMenu $f.msg -label "Print Messages : " \
	    -options {
		label.width  19
		label.anchor e
		menubutton.width 15
	    } 

	foreach step {0 1 2 3 4 5 } {
	    $f.msg add command $step -label $msg($step)
	}

	$f.msg config -variable options.printmsg 
	pack $f.msg
	





	

	checkbutton $f.parthread -text "Parallel meshing thread" \
	    -variable options.parthread
	checkbutton $f.second -text "Second order elements" \
	    -variable options.secondorder
	checkbutton $f.quad -text "Quad dominated" \
	    -variable options.quad -command {
		if { ${options.quad} } {
		    set meshoptions.laststep os
		}
	    }
	checkbutton $f.invtets -text "Invert volume elements" \
	    -variable options.inverttets
	checkbutton $f.invtrigs -text "Invert surface elements" \
	    -variable options.inverttrigs
	checkbutton $f.azref -text "Automatic Z-refinement" \
	    -variable options.autozrefine

	pack $f.parthread $f.second $f.quad $f.invtets $f.invtrigs $f.azref 



	tixControl $f.elementorder -label "Element order: " -integer true \
	    -variable options.elementorder -min 1 -max 20 \
	    -options {
		entry.width 2
		label.width 20
		label.anchor e
	    }	

        pack $f.elementorder


#	tixControl $f.memory -label "Large Memory \[MB\]: " -integer true \
\#	    -variable options.memory -min 0 -max 2000 \
\#	    -options {
#		entry.width 5
#		label.width 20
#		label.anchor e
#	    }	

#	global userlevel
#	if { $userlevel >= 3} { pack $f.memory }


	# Mesh - Size options
	set f [$w.nb subwidget meshsize]

	tixControl $f.meshsize -label "max mesh-size: " -integer false \
	    -variable options.meshsize -min 1e-9 -max 1e6 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	

	tixControl $f.minmeshsize -label "min mesh-size: " -integer false \
	    -variable options.minmeshsize -min 0 -max 1e6 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	

	tixControl $f.grading -label "mesh-size grading: " -integer false \
	    -variable options.grading -min 0.1 -max 1 -step 0.1 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	
	
	pack $f.meshsize $f.minmeshsize $f.grading



	frame $f.msf 
	pack $f.msf

 	tixLabelEntry $f.msf.ent -label "mesh-size file: "  \
 	    -labelside top \
 	    -options {  
 		entry.textVariable options.meshsizefilename 
 		entry.width 35
 		label.width 25
 		label.anchor w
 	    }	
 	button $f.msf.btn -text "Browse" -command {
	    global options.meshsizefilename
	    set types {
		{"Meshsize file"   {.msz}	} }
	    set options.meshsizefilename [tk_getOpenFile -filetypes $types -initialfile ${options.meshsizefilename}]
	}

 	pack $f.msf.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $f.msf.btn -side left -anchor s -padx 4 -pady 4


	label $f.lab -text "Additional mesh size restrictions:"

	#csg-meshsize options

	frame $f.csg -relief groove -borderwidth 3
	pack $f.csg -fill x


	frame $f.csg.curv
	pack $f.csg.curv -anchor w

	scale $f.csg.curv.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable options.curvaturesafety
	label $f.csg.curv.la -text "Elements per curvature radius"
	pack $f.csg.curv.sc $f.csg.curv.la -side left 

	frame $f.csg.elen
	pack $f.csg.elen -anchor w
	scale $f.csg.elen.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable options.segmentsperedge
	label $f.csg.elen.la -text "Elements per edge"
	pack $f.csg.elen.sc $f.csg.elen.la -side left
	

	#stl-meshsize options

	frame $f.stl -relief groove -borderwidth 3
	pack $f.stl -fill x

	frame $f.stl.r2
	pack $f.stl.r2 -anchor w
	scale $f.stl.r2.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable stloptions.resthchartdistfac
	checkbutton $f.stl.r2.bu -text "STL - chart distance" \
	    -variable stloptions.resthchartdistenable
	pack $f.stl.r2.sc $f.stl.r2.bu -side left
	
	frame $f.stl.r6
	pack $f.stl.r6 -anchor w
	scale $f.stl.r6.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable stloptions.resthlinelengthfac
	checkbutton $f.stl.r6.bu -text "STL - line length" \
	    -variable stloptions.resthlinelengthenable
	pack $f.stl.r6.sc $f.stl.r6.bu -side left
	
	frame $f.stl.r3
	pack $f.stl.r3 -anchor w
	scale $f.stl.r3.sc -orient horizontal -length 200 -from 0.2 -to 8 \
	    -resolution 0.1 -variable stloptions.resthcloseedgefac 
	checkbutton $f.stl.r3.bu -text "STL/IGES/STEP - close edges" \
	    -variable stloptions.resthcloseedgeenable 
	
	pack $f.stl.r3.sc $f.stl.r3.bu -side left
	
	frame $f.stl.r1
	pack $f.stl.r1 -anchor w
	scale $f.stl.r1.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable stloptions.resthsurfcurvfac
	checkbutton $f.stl.r1.bu -text "STL - surface curvature" \
	    -variable stloptions.resthsurfcurvenable
	pack $f.stl.r1.sc $f.stl.r1.bu -side left

	frame $f.stl.r3b
	pack $f.stl.r3b -anchor w
	scale $f.stl.r3b.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable stloptions.resthedgeanglefac
	checkbutton $f.stl.r3b.bu -text "STL - edge angle" \
	-variable stloptions.resthedgeangleenable
	pack $f.stl.r3b.sc $f.stl.r3b.bu -side left
	
	frame $f.stl.r5
	pack $f.stl.r5 -anchor w
	scale $f.stl.r5.sc -orient horizontal -length 200 -from 0.2 -to 5 \
	    -resolution 0.1 -variable stloptions.resthsurfmeshcurvfac
	checkbutton $f.stl.r5.bu -text "STL - surface mesh curv" \
	    -variable stloptions.resthsurfmeshcurvenable
	pack $f.stl.r5.sc $f.stl.r5.bu -side left
	
	
	checkbutton $f.stl.recalch -text "STL - Recalc mesh size for surface optimization" \
	    -variable stloptions.recalchopt
	pack $f.stl.recalch

	button $f.stl.calch -text "Calc New H" -command { redraw; Ng_STLCalcLocalH }
	pack $f.stl.calch






	set f [$w.nb subwidget chartopt]


	label $f.lab1 -text "Yellow Edges Angle ()"
	scale $f.scale1 -orient horizontal -length 300 \
	    -from 0 -to 90 -resolution 1  -tickinterval 10 \
	    -variable  stloptions.yangle 

	pack $f.lab1 $f.scale1

	label $f.lab2e -text "Edge Corner Angle ()"
	scale $f.scale2e -orient horizontal -length 360 -from 0 -to 180 \
	    -resolution 1  -tickinterval 20 \
	    -variable  stloptions.edgecornerangle 
	pack $f.lab2e $f.scale2e
	
	label $f.lab2 -text "Chart Angle ()"
	scale $f.scale2 -orient horizontal -length 360 -from 0 -to 180 \
	    -resolution 1  -tickinterval 20 \
	    -variable  stloptions.chartangle 
	pack $f.lab2 $f.scale2
	
	label $f.lab2b -text "Outer Chart Angle ()"
	scale $f.scale2b -orient horizontal -length 360 -from 0 -to 180 \
	    -resolution 1  -tickinterval 20 \
	    -variable  stloptions.outerchartangle 
	pack $f.lab2b $f.scale2b




	# Optimization options
	
	set f [$w.nb subwidget optimizer]

	


	tixControl $f.os2d -label "Surface opt steps: " -integer true \
	    -variable options.optsteps2d -min 0 -max 99 -step 1 \
	    -options {
		entry.width 3
		label.width 25
		label.anchor e
	    }	

	tixControl $f.os3d -label "Volume opt steps: " -integer true \
	    -variable options.optsteps3d -min 0 -max 99 -step 1 \
	    -options {
		entry.width 3
		label.width 25
		label.anchor e
	    }	
	
	tixControl $f.elw -label "Element size weight: " -integer false \
	    -variable options.elsizeweight -min 0 -max 1 -step 0.1 \
	    -options {
		entry.width 3
		label.width 25
		label.anchor e
	    }	

	tixControl $f.wem -label "Worst element measure: " -integer false \
	    -variable options.opterrpow -min 1 -max 10 -step 1 \
	    -options {
		entry.width 3
		label.width 25
		label.anchor e
	    }	

	pack $f.os2d $f.os3d $f.elw $f.wem
	
	frame $f.badellimit
	pack $f.badellimit -fill x
	label $f.badellimit.lab -text "bad element criterion";
	scale $f.badellimit.scale -orient horizontal -length 150 \
	    -from 160 -to 180 -resolution 1 \
	-variable options.badellimit
	pack $f.badellimit.scale $f.badellimit.lab -side right -anchor s


	# insider options
	set f [$w.nb subwidget insider]
	


	checkbutton $f.localh -text "Use Local Meshsize" \
	    -variable options.localh
	checkbutton $f.delauney -text "Use Delaunay" \
	    -variable options.delaunay
	checkbutton $f.checkoverlap -text "Check Overlapping" \
	    -variable options.checkoverlap
	checkbutton $f.checkcb -text "Check Chart Boundary" \
	    -variable options.checkchartboundary
	checkbutton $f.blockfill -text "Do Blockfilling" \
	    -variable options.blockfill

	pack  $f.localh  $f.delauney $f.checkoverlap  $f.blockfill $f.checkcb   -anchor w



	
	# debugging options
	set f [$w.nb subwidget debug]

	frame $f.cb
	pack $f.cb -side top

	

	checkbutton $f.cb.slowchecks -text "Slow checks" \
	    -variable debug.slowchecks -command { Ng_SetDebugParameters }
	checkbutton $f.cb.debugoutput -text "Debugging outout" \
	    -variable debug.debugoutput -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltexline -text "Halt on exising line" \
	    -variable debug.haltexistingline  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltoverlap -text "Halt on Overlap" \
	    -variable debug.haltoverlap  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltsuc -text "Halt on success" \
	    -variable debug.haltsuccess  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltnosuc -text "Halt on no success" \
	    -variable debug.haltnosuccess  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltlargequal -text "Halt on large quality class" \
	    -variable debug.haltlargequalclass  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltseg -text "Halt on Segment:" \
	    -variable debug.haltsegment  -command { Ng_SetDebugParameters }
	checkbutton $f.cb.haltnode -text "Halt on Node:" \
	    -variable debug.haltnode  -command { Ng_SetDebugParameters }


	pack $f.cb.slowchecks $f.cb.debugoutput $f.cb.haltexline $f.cb.haltoverlap $f.cb.haltsuc $f.cb.haltnosuc $f.cb.haltlargequal  $f.cb.haltseg   $f.cb.haltnode 

	frame $f.cb.hf
	pack $f.cb.hf -pady 5
	checkbutton $f.cb.hf.cb -text "Halt on Face:" \
	    -variable debug.haltface  -command { Ng_SetDebugParameters }
	entry $f.cb.hf.ent -textvariable debug.haltfacenr -width 5  
	pack $f.cb.hf.cb $f.cb.hf.ent -side left 

	checkbutton $f.cb.showactivechart -text "Show Active Meshing-Chart" \
	    -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }




	pack $f.cb.showactivechart
	

	frame $f.segs
	pack $f.segs -pady 5
	label $f.segs.lab1 -text "P1:";
	entry $f.segs.ent1 -width 8 -relief sunken \
	    -textvariable debug.haltsegmentp1 
	label $f.segs.lab2 -text "P2:";
	entry $f.segs.ent2 -width 8 -relief sunken \
	    -textvariable debug.haltsegmentp2 
	pack $f.segs.lab1 $f.segs.ent1 $f.segs.lab2 $f.segs.ent2  -side left



	frame $f.cont -relief groove -borderwidth 3
	pack $f.cont 
	#-fill x 
	
	checkbutton $f.cont.multidrawing -text "Draw Meshing" \
	    -variable multithread_drawing 
	pack $f.cont.multidrawing
	
	checkbutton $f.cont.multitestmode -text "Meshing Testmode" \
	-variable multithread_testmode 
	pack $f.cont.multitestmode
	
	button $f.cont.goon -text "Go On" -command { set multithread_pause 0 }
	pack $f.cont.multidrawing $f.cont.multitestmode $f.cont.goon -side left -expand yes
	



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
	
# 	pack $w.bbox -side bottom -fill x

	
	frame $w.bu
	pack $w.bu -fill x -ipady 3

	button $w.bu.apl -text "Apply" -command { 
	    [.options_dlg.nb subwidget meshsize].meshsize invoke
	    [.options_dlg.nb subwidget meshsize].grading invoke
	    [.options_dlg.nb subwidget optimizer].os2d invoke
	    [.options_dlg.nb subwidget optimizer].os3d invoke
	    [.options_dlg.nb subwidget optimizer].elw invoke
	    [.options_dlg.nb subwidget optimizer].wem invoke

	    Ng_SetMeshingParameters 
	    Ng_SetDebugParameters
	}

	button $w.bu.ok -text "Done" -command {
	    [.options_dlg.nb subwidget meshsize].meshsize invoke
	    [.options_dlg.nb subwidget meshsize].grading invoke
	    [.options_dlg.nb subwidget optimizer].os2d invoke
	    [.options_dlg.nb subwidget optimizer].os3d invoke
	    [.options_dlg.nb subwidget optimizer].elw invoke
	    [.options_dlg.nb subwidget optimizer].wem invoke

	    Ng_SetMeshingParameters
	    Ng_SetDebugParameters
	    wm withdraw .options_dlg
#	    destroy .options_dlg
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
 



	tixNoteBook $w.nb -ipadx 6 -ipady 6
	
	$w.nb add general -label "General" -underline 0
	$w.nb add stl -label "STL" -underline 0
	$w.nb add occ -label "IGES/STEP" -underline 0
	$w.nb add mesh -label "Mesh"   -underline 0
	$w.nb add light -label "Light"   -underline 0
	$w.nb add edges -label "Edges"   -underline 0
	$w.nb add misc -label "Misc."  -underline 3


	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	





	# general
	set f [$w.nb subwidget general]

	checkbutton $f.backcol -text "White Background" \
	-variable viewoptions.whitebackground \
	-command { Ng_SetVisParameters; redraw }

	checkbutton $f.cross -text "Draw Coordinate Cross" \
	-variable viewoptions.drawcoordinatecross \
	-command { Ng_SetVisParameters; redraw }

	checkbutton $f.color -text "Draw Color-bar" \
	-variable viewoptions.drawcolorbar \
	-command { Ng_SetVisParameters; redraw }

	checkbutton $f.netgen -text "Draw Netgen-logo" \
	-variable viewoptions.drawnetgenlogo \
	-command { Ng_SetVisParameters; redraw }

	pack $f.backcol $f.cross $f.color $f.netgen

# 	checkbutton $f.stereo -text "Stereo View" \
# 	-variable viewoptions.stereo \
# 	-command { Ng_SetVisParameters; redraw }
# 	pack $f.stereo

	# stl geometry 
	set f [$w.nb subwidget stl]

	frame $f.show -relief groove -borderwidth 3
	pack $f.show
	checkbutton $f.show.showtrias -text "Show STL-Triangles" \
	    -variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }
	pack $f.show.showtrias -anchor w
	
	checkbutton $f.show.showfilledtrias -text "Show Filled Triangles" \
	    -variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }
	pack $f.show.showfilledtrias -anchor w
	
	checkbutton $f.show.showactivechart -text "Show Active Meshing-Chart" \
	    -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }
	pack $f.show.showactivechart -anchor w
	
	checkbutton $f.show.showedges -text "Show Edges" \
	    -variable stloptions.showedges -command { Ng_SetVisParameters; redraw }
	pack $f.show.showedges -anchor w
	
	frame $f.special -relief groove -borderwidth 3
	pack $f.special
	checkbutton $f.special.showmarktrias -text "Show Chart Triangles" \
	    -variable stloptions.showmarktrias \
	    -command {set stldoctor.showfaces 0; Ng_STLDoctor; Ng_SetVisParameters; redraw }
	pack $f.special.showmarktrias -side left

	checkbutton $f.special.showfaces -text "Show Faces" \
	    -variable stldoctor.showfaces \
	    -command {set stloptions.showmarktrias 0; Ng_STLDoctor; Ng_SetVisParameters; redraw}    
	pack $f.special.showfaces -side left

	frame $f.fn -relief groove -borderwidth 3
	pack $f.fn
	label $f.fn.lab3 -text "Chart/Face number:"
	scale $f.fn.scale3 -orient horizontal -length 200 -from 0 -to 200 \
	    -resolution 1  -tickinterval 50 \
	    -command { Ng_SetVisParameters; redraw } -variable  stloptions.chartnumber 
	pack $f.fn.lab3 $f.fn.scale3 -side left
	
	frame $f.fo -relief groove -borderwidth 3
	pack $f.fo
	label $f.fo.lab -text "Chart/Face Offset:";
	entry $f.fo.ent -width 5 -relief sunken \
	    -textvariable stloptions.chartnumberoffset
	pack $f.fo.lab $f.fo.ent -side left

	frame $f.mt
	pack $f.mt -fill x
	checkbutton $f.mt.bu -text "Show Marked (Dirty) Triangles" \
	    -variable stldoctor.showmarkedtrigs \
	    -command {Ng_STLDoctor; redraw}    
	pack $f.mt.bu

	frame $f.ep
	pack $f.ep -fill x
	checkbutton $f.ep.bu -text "show edge corner points" \
	    -variable stldoctor.showedgecornerpoints \
	    -command {Ng_STLDoctor; redraw}    
	pack $f.ep.bu

	frame $f.stt
	pack $f.stt -fill x
	checkbutton $f.stt.bu -text "show touched triangle chart" \
	    -variable stldoctor.showtouchedtrigchart \
	    -command {set stldoctor.showfaces 0; set stloptions.showmarktrias 1; \
			  Ng_STLDoctor; Ng_SetVisParameters; redraw}    
	pack $f.stt.bu

	frame $f.sml
	pack $f.sml -fill x
	checkbutton $f.sml.bu -text "draw meshed edges" \
	    -variable stldoctor.drawmeshededges \
	    -command {Ng_STLDoctor;}    
	pack $f.sml.bu
	
	
	frame $f.sm
	pack $f.sm -fill x
	checkbutton $f.sm.bu -text "select with mouse" \
	    -variable stldoctor.selectwithmouse
	pack $f.sm.bu
	
	frame $f.st -relief groove -borderwidth 3
	pack $f.st -fill x
	label $f.st.lab -text "Select triangle by number";
	entry $f.st.ent -width 5 -relief sunken \
	    -textvariable stldoctor.selecttrig
	pack $f.st.ent $f.st.lab -side left -expand yes
	
	frame $f.vc -relief groove -borderwidth 3
	pack $f.vc -fill x
	checkbutton $f.vc.bu -text "show vicinity" \
	    -variable stldoctor.showvicinity \
	    -command {Ng_STLDoctor vicinity; redraw}
	label $f.vc.lab -text "vicinity size";
	scale $f.vc.sc -orient horizontal -length 200 -from 0 -to 200 \
	    -resolution 1 -variable stldoctor.vicinity \
	    -command { Ng_STLDoctor vicinity; redraw }
	pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes
	


	# IGES/STEP
	set f [$w.nb subwidget occ]
	
	checkbutton $f.occshowsurfaces -text "Show surfaces " \
	    -variable occoptions.showsurfaces \
	    -command { Ng_SetOCCVisParameters; redraw }

	checkbutton $f.occshowedges -text "Show edges " \
	    -variable occoptions.showedges \
	    -command { Ng_SetOCCVisParameters; redraw }

	frame $f.deflection -relief groove -borderwidth 3
	pack $f.deflection -fill x
	button $f.deflection.lab -text "Rebuild visualization data" \
	    -command {
		Ng_SetOCCVisParameters
		Ng_OCCCommand buildvisualizationmesh
		redraw
	    }

	tixControl $f.deflection.ent -label "Visualization smoothness" -integer false \
	    -variable occoptions.deflection -min 0.1 -max 3 -step 0.1 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters }





	pack $f.deflection.ent $f.deflection.lab -side left  -expand yes
	pack $f.occshowsurfaces $f.occshowedges


	# ACIS visualization / construction

	tixControl $f.showsolid -label "Show solid (0 for all)" -integer true \
            -variable occoptions.showsolidnr -min 0 -max 999 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters; redraw }
    
	tixControl $f.showsolid2 -label "Show solid 2" -integer true \
            -variable occoptions.showsolidnr2 -min 0 -max 999 \
	    -options { entry.width 3 } \
	    -command { Ng_SetOCCVisParameters; redraw }

	button $f.subtract -text "Subtract (2 minus 1)" \
	    -command {
		Ng_ACISCommand subtract ${occoptions.showsolidnr} ${occoptions.showsolidnr2}
		redraw
	    }

	button $f.combine -text "Combine all" \
	    -command {
		Ng_ACISCommand combineall
		redraw
	    }

	pack $f.showsolid $f.showsolid2 $f.subtract $f.combine




	# mesh options
	set f [$w.nb subwidget mesh]

	checkbutton $f.showcolor -text "Colored Meshsize Visualization" \
	    -variable viewoptions.colormeshsize \
	    -command { Ng_SetVisParameters; redraw }


	checkbutton $f.showfilledtrigs -text "Show filled triangles" \
	-variable viewoptions.drawfilledtrigs \
	-command { Ng_SetVisParameters; redraw }
	
	checkbutton $f.showedges -text "Show edges" \
	-variable viewoptions.drawedges \
	-command { Ng_SetVisParameters; redraw }
	

	checkbutton $f.showoutline -text "Show Triangle Outline" \
	    -variable viewoptions.drawoutline \
	    -command { Ng_SetVisParameters; redraw }

	
	tixControl $f.subdiv -label "Subdivision" -integer true \
            -variable visoptions.subdivisions -min 0 -max 8 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; Ng_Vis_Set parameters; Ng_SetNextTimeStamp; redraw }
	
	
	checkbutton $f.showbadels -text "Show bad elements" \
	    -variable viewoptions.drawbadels \
	    -command { Ng_SetVisParameters; redraw }



	checkbutton $f.showprisms -text "Show prisms" \
	    -variable viewoptions.drawprisms \
	    -command { Ng_SetVisParameters; redraw }

	checkbutton $f.showpyramids -text "Show pyramids" \
	    -variable viewoptions.drawpyramids \
	    -command { Ng_SetVisParameters; redraw }

	checkbutton $f.showhexes -text "Show hexes" \
	    -variable viewoptions.drawhexes \
	    -command { Ng_SetVisParameters; redraw }

	frame $f.fshrink
	label $f.fshrink.lab -text "Shrink elements"
	scale $f.fshrink.scale -orient horizontal -length 200 -from 0 -to 1.0001 \
	    -resolution 0.01  -tickinterval 0.25 \
	    -command { Ng_SetVisParameters; after idle redraw } \
            -variable  viewoptions.shrink


	checkbutton $f.showidentified -text "Show identified points" \
	    -variable viewoptions.drawidentified \
	    -command { Ng_SetVisParameters; redraw }

	checkbutton $f.showmetispartition -text "Show METIS Partition" \
	    -variable viewoptions.drawmetispartition \
	    -command { Ng_SetVisParameters; redraw }

	checkbutton $f.showpointnumbers -text "Show Point-numbers" \
	    -variable viewoptions.drawpointnumbers \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showedgenumbers -text "Show Edge-numbers" \
	    -variable viewoptions.drawedgenumbers \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showfacenumbers -text "Show Face-numbers" \
	    -variable viewoptions.drawfacenumbers \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showelementnumbers -text "Show Element-numbers" \
	    -variable viewoptions.drawelementnumbers \
	    -command { Ng_SetVisParameters; redraw }
	
	# label $f.showdomainlab -text "Domain Surface"
#	scale $f.showdomain -orient horizontal -length 100 -from 0 -to 50 \
	    -resolution 1 -variable  viewoptions.drawdomainsurf    \
	    -command { Ng_SetVisParameters; redraw } \
	    -label "Domain Surface" 

	tixControl $f.showdomain -label "Show surface of domain" -integer true \
            -variable viewoptions.drawdomainsurf -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; Ng_Vis_Set parameters; redraw }
    
    
    

	frame $f.center -relief groove -borderwidth 3
	pack $f.center -fill x
	button $f.center.lab -text "Set Center Point" \
	    -command { Ng_SetVisParameters; Ng_Center; redraw }
	entry $f.center.ent -width 5 -relief sunken \
	    -textvariable viewoptions.centerpoint 
	pack $f.center.ent $f.center.lab -side left  -expand yes
	
	frame $f.drawel -relief groove -borderwidth 3
	pack $f.drawel -fill x
	button $f.drawel.lab -text "Draw Element" \
	    -command { Ng_SetVisParameters; Ng_ZoomAll; redraw }
	entry $f.drawel.ent -width 5 -relief sunken \
	    -textvariable viewoptions.drawelement 
	pack $f.drawel.ent $f.drawel.lab -side left  -expand yes

        pack $f.showfilledtrigs
	pack $f.showoutline $f.subdiv $f.showedges  $f.showbadels 
	# pack $f.showdomainlab 
	pack $f.showdomain 
	pack $f.showpointnumbers 
	pack $f.showedgenumbers $f.showfacenumbers $f.showelementnumbers 
	pack $f.showmetispartition


	frame $f.frametets
	checkbutton $f.frametets.showtets -text "Show Tets in domain " \
	    -variable viewoptions.drawtets \
	    -command { Ng_SetVisParameters; redraw }
	tixControl $f.frametets.showtetsdomain -label "" -integer true \
	    -variable viewoptions.drawtetsdomain -min 0 -max 500 \
	    -options { entry.width 2 } \
	    -command { Ng_SetVisParameters; redraw }

	pack $f.frametets
	pack $f.frametets.showtets $f.frametets.showtetsdomain -side left


	pack $f.showcolor    $f.showpyramids $f.showprisms $f.showhexes $f.showidentified
	
	pack $f.fshrink 
	pack $f.fshrink.lab $f.fshrink.scale -side left
	
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
	set f [$w.nb subwidget light]
	
	label $f.lab1 -text "Ambient Light"
	scale $f.scale1 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.amb 
	label $f.lab2 -text "Diffuse Light"
	scale $f.scale2 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.diff 
	label $f.lab3 -text "Specular Light"
	scale $f.scale3 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.spec 
	label $f.lab4 -text "Material Shininess"
	scale $f.scale4 -orient horizontal -length 300 -from 0 -to 128 \
	    -resolution 1  -tickinterval 32 \
	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.mat.shininess 
	label $f.lab5 -text "Material Transparency"
	scale $f.scale5 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	-command { Ng_SetVisParameters; redraw } -variable  viewoptions.mat.transp 
	
	pack $f.lab1 $f.scale1 $f.lab2 $f.scale2 $f.lab3 $f.scale3 $f.lab4 $f.scale4 $f.lab5 $f.scale5
	




	# edges options
	set f [$w.nb subwidget edges]

	checkbutton $f.showedges -text "Show Edges" \
	    -variable viewoptions.drawededges \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showpoints -text "Show Points" \
	    -variable viewoptions.drawedpoints \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showpointnrs -text "Show Points Nrs" \
	    -variable viewoptions.drawedpointnrs \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.showtang -text "Show CP Tangents" \
	    -variable viewoptions.drawedtangents \
	    -command { Ng_SetVisParameters; redraw }
	checkbutton $f.drawedgenrs -text "Show Edge Nrs" \
	    -variable viewoptions.drawededgenrs \
	    -command { Ng_SetVisParameters; redraw }
	
	pack $f.showedges $f.showpoints $f.showpointnrs $f.showtang $f.drawedgenrs

	frame $f.center -relief groove -borderwidth 3
	pack $f.center -fill x
	button $f.center.lab -text "Set Center Point" \
	    -command { Ng_SetVisParameters; Ng_Center; redraw }
	entry $f.center.ent -width 5 -relief sunken \
	    -textvariable viewoptions.centerpoint 
	pack $f.center.ent $f.center.lab -side left  -expand yes
	


	frame $f.f1
	pack $f.f1 -pady 5
	label $f.f1.lab -text "SpecPoint Veclen"
	entry $f.f1.ent -width 5 -relief sunken -textvariable viewoptions.specpointvlen
	pack $f.f1.lab $f.f1.ent
	



	# misc options
	set f [$w.nb subwidget misc]

	frame $f.point -relief groove -borderwidth 3

	frame $f.point.dp
	
	checkbutton $f.point.dp.drawpoint -text "Draw Point" \
	    -variable viewoptions.drawspecpoint \
	    -command { Ng_SetVisParameters; redraw }

	entry $f.point.dp.px -width 8 -relief sunken -textvariable viewoptions.specpointx
	entry $f.point.dp.py -width 8 -relief sunken -textvariable viewoptions.specpointy
	entry $f.point.dp.pz -width 8 -relief sunken -textvariable viewoptions.specpointz

	pack $f.point.dp.drawpoint $f.point.dp.px $f.point.dp.py $f.point.dp.pz -side left

	pack $f.point.dp

	checkbutton $f.point.center -text "Use as Center" \
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

	pack $f.point.center
	
	pack $f.point -fill x -ipady 3


	
	frame $w.bu
	pack $w.bu -fill x -ipady 3


	button $w.bu.done -text "Done" -command {
	    Ng_SetVisParameters;
	    redraw
	    destroy .viewopts_dlg
	}
	button $w.bu.apply -text "Apply" -command {
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

	label $w.lab1 -text "Normal x"
	scale $w.scale1 -orient horizontal -length 300 -from -1 -to 1 \
	    -resolution 0.01  -tickinterval 0.5 \
	    -variable  viewoptions.clipping.nx \
	    -command { clipplanecommand }
#	    -command { popupcheckredraw2 clippingdialog_pop1 ${viewoptions.clipping.enable} }

#		Ng_SetVisParameters; 
#		if { ${viewoptions.clipping.enable} == 1 } { redraw };
#		Ng_SetVisParameters 
	
	label $w.lab2 -text "Normal y"
	scale $w.scale2 -orient horizontal -length 300 -from -1 -to 1 \
	    -resolution 0.01  -tickinterval 0.5 \
	    -variable  viewoptions.clipping.ny \
	    -command { clipplanecommand }
#	    -command { popupcheckredraw2 clippingdialog_pop2 ${viewoptions.clipping.enable} }

	label $w.lab3 -text "Normal z"
	scale $w.scale3 -orient horizontal -length 300 -from -1 -to 1 \
	    -resolution 0.01  -tickinterval 0.5 \
	    -variable  viewoptions.clipping.nz \
	    -command { clipplanecommand }
#	    -command { popupcheckredraw2 clippingdialog_pop3 ${viewoptions.clipping.enable} }
	label $w.lab4 -text "Distance"
	scale $w.scale4 -orient horizontal -length 300 -from -1 -to 1.001 \
	    -resolution 0.0001  -tickinterval 0.5 \
	    -variable  viewoptions.clipping.dist \
	    -command { clipplanecommand }
#	    -command { popupcheckredraw2 clippingdialog_pop4 ${viewoptions.clipping.enable} }

	label $w.lab5 -text "Additional Distance"
	scale $w.scale5 -orient horizontal -length 300 -from -1 -to 1.001 \
	    -resolution 0.0001  -tickinterval 0.5 \
	    -variable  viewoptions.clipping.dist2 \
	    -command { clipplanecommand }

	
	
	tixControl $w.clipdomain -label "Clip only domain" -integer true \
	    -variable viewoptions.clipping.onlydomain -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { clipplanecommand; }
#	    -command { Ng_SetVisParameters; redraw }
	tixControl $w.donotclipdomain -label "Do not clip domain" -integer true \
	    -variable viewoptions.clipping.notdomain -min 0 -max 50 \
	    -options { entry.width 2 } \
	    -command { clipplanecommand; }
#	    -command { Ng_SetVisParameters; redraw }

	pack $w.lab1 $w.scale1 $w.lab2 $w.scale2 $w.lab3 $w.scale3 $w.lab4 $w.scale4 $w.lab5 $w.scale5 $w.clipdomain $w.donotclipdomain

	
	checkbutton $w.cb1 -text "Enable clipping" \
	    -variable viewoptions.clipping.enable \
	    -command { Ng_SetVisParameters; redraw } 
	
	pack $w.cb1
	

	
	frame $w.bu
#	pack $w.bu -fill x
	pack $w.bu -fill x -ipady 3

	button $w.bu.cancle -text "Done" -command "destroy $w"
	pack $w.bu.cancle  -expand yes
	
	
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

	
	tixControl $w.meshsize -label "max mesh-size: " -integer false \
	    -variable options.meshsize -min 1e-6 -max 1e6 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	

	pack $w.meshsize

	global localh
	set localh 1
	tixControl $w.loch -label "local mesh-size: " -integer false \
	    -variable localh -min 1e-6 -max 1e6 \
	    -options {
		entry.width 6
		label.width 25
		label.anchor e
	    }	
	
	pack $w.loch
	
	
	button $w.restface -text "Restrict H at face"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH face $localh
	    }
	button $w.restedge -text "Restrict H at edge"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH edge $localh
	    }
	button $w.restelement -text "Restrict H at element"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH element $localh
	    }
	button $w.restpoint -text "Restrict H at point"  \
	    -command {
		.refinement_dlg.meshsize invoke
		.refinement_dlg.loch invoke
		Ng_RestrictH point $localh
	    }


	pack $w.restface $w.restedge $w.restelement $w.restpoint



	button $w.anisoedge -text "Declare Anisotropic edge"  \
	    -command {
		Ng_Anisotropy edge 
	    }
	pack $w.anisoedge
	

	frame $w.bu
	pack $w.bu -fill x -ipady 3


	button $w.bu.cancle -text "Done" -command "destroy .refinement_dlg"
	button $w.bu.refine -text "Refine"  \
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
	button $w.bu.zrefine -text "Z-Refine"  \
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

	tixNoteBook $w.nb -ipadx 6 -ipady 6
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

    set wd .stldoctor_dlg

    if {[winfo exists .stldoctor_dlg] == 1} {
	wm withdraw $wd
	wm deiconify $wd
	focus $wd 
    } {
	
    toplevel $wd

    tixNoteBook $wd.nb -ipadx 6 -ipady 6
	
    $wd.nb add general -label "General" -underline 0
    $wd.nb add topology -label "Edit Topology"  -underline 5
    $wd.nb add edges -label "Edit Edges"   -underline 5
    $wd.nb add normals -label "Edit Normals"   -underline 5
    $wd.nb add advanced -label "Advanced"   -underline 0


    pack $wd.nb -expand yes -fill both -padx 5 -pady 5 -side top	


    # GENERAL *****************************

    set f [$wd.nb subwidget general]


    frame $f.show
    pack $f.show -fill x
    checkbutton $f.show.showtrias -text "Show STL-Triangles" \
	-variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }
    pack $f.show.showtrias -anchor w
    
    checkbutton $f.show.showfilledtrias -text "Show Filled Triangles" \
	-variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }
    pack $f.show.showfilledtrias -anchor w

    set selmodevals { 0 1 2 3 4 }
    set selmodelabs(0) "triangle" 
    set selmodelabs(1) "edge" 
    set selmodelabs(2) "point" 
    set selmodelabs(3) "line" 
    set selmodelabs(4) "line cluster" 

    tixOptionMenu $f.selmode -label "Double Click selects :" \
	-options {
	    label.width  19
	    label.anchor e
	    menubutton.width 15
	} 

    foreach selmodev $selmodevals {
	$f.selmode add command $selmodev -label $selmodelabs($selmodev)
    }
    $f.selmode config -variable stldoctor.selectmode
    $f.selmode config -command { Ng_STLDoctor }
    global stldoctor.selectmode
    pack $f.selmode

    frame $f.sm
    pack $f.sm -fill x
    checkbutton $f.sm.bu -text "select with mouse" \
	-variable stldoctor.selectwithmouse
    pack $f.sm.bu 

    frame $f.st -relief groove -borderwidth 3
    pack $f.st -fill x
    label $f.st.lab -text "Select triangle by number";
    entry $f.st.ent -width 5 -relief sunken \
	-textvariable stldoctor.selecttrig
    pack $f.st.ent $f.st.lab -side left -expand yes

    frame $f.vc -relief groove -borderwidth 3
    pack $f.vc -fill x
    checkbutton $f.vc.bu -text "show vicinity" \
	-variable stldoctor.showvicinity \
	-command {Ng_STLDoctor vicinity; redraw}
    label $f.vc.lab -text "vicinity size";
    scale $f.vc.sc -orient horizontal -length 200 -from 0 -to 200 \
	-resolution 1 -variable stldoctor.vicinity \
	-command { Ng_STLDoctor vicinity; redraw }
    pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes

    frame $f.ge -relief groove -borderwidth 3
    pack $f.ge -fill x
    button $f.ge.neighbourangles -text "calc neighbourangles" -command {Ng_STLDoctor neighbourangles}
    button $f.ge.showcoords -text "show coords of touched triangle" -command {Ng_STLDoctor showcoords}
    button $f.ge.moveptm -text "move point to middle of trianglepoints" -command {Ng_STLDoctor movepointtomiddle; redraw}
    button $f.ge.destroy0trigs -text "destroy 0-volume triangles" -command {Ng_STLDoctor destroy0trigs}
    pack $f.ge.neighbourangles $f.ge.showcoords $f.ge.moveptm $f.ge.destroy0trigs -expand yes 


    button $f.ge.cancle -text "Done" -command {destroy .stldoctor_dlg }
    pack $f.ge.cancle -expand yes

    # TOPOLOGY ********************
    set f [$wd.nb subwidget topology]

    frame $f.oc -relief groove -borderwidth 3
    pack $f.oc -fill x
    button $f.oc.bu -text "invert orientation of selected trig" -command {Ng_STLDoctor invertselectedtrig; redraw }
    button $f.oc.bu2 -text "orient after selected trig" -command {Ng_STLDoctor orientafterselectedtrig; redraw }
    pack $f.oc.bu $f.oc.bu2 -side left  -expand yes

    button $f.toperr -text "mark inconsistent triangles" -command {Ng_STLDoctor marktoperrortrigs; redraw }

    button $f.deltrig -text "delete selected triangle" -command {Ng_STLDoctor deleteselectedtrig; redraw }
    button $f.geosmooth -text "geometric smoothing" -command {Ng_STLDoctor smoothgeometry; redraw }

    pack $f.toperr $f.deltrig $f.geosmooth





    # EDGES ***********************
    set f [$wd.nb subwidget edges]


    frame $f.be -relief groove -borderwidth 3 
    pack $f.be -fill x
    label $f.be.lab -text "build edges with yellow angle:";
    scale $f.be.sc -orient horizontal -length 200 -from 0 -to 100 \
	-resolution 0.5
    $f.be.sc config -variable stloptions.yangle 
    $f.be.sc config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }
    label $f.be.lab2 -text "continue edges with yellow angle:";
    scale $f.be.sc2 -orient horizontal -length 200 -from 0 -to 100 \
	-resolution 0.5
    $f.be.sc2 config -variable stloptions.contyangle 
    $f.be.sc2 config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }



    button $f.be.buildedges -text "Build Edges" -command {Ng_STLDoctor buildedges; redraw}
    pack $f.be.lab $f.be.sc $f.be.lab2 $f.be.sc2 $f.be.buildedges -expand yes

    frame $f.se
    pack $f.se -fill x
    checkbutton $f.se.bu -text "show excluded" \
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

    tixOptionMenu $f.edgeselmode -label "Double Click sets edge :" \
	-options {
	    label.width  19
	    label.anchor e
	    menubutton.width 15
	} 

    foreach edgeselmodev $edgeselmodevals {
	$f.edgeselmode add command $edgeselmodev -label $edgeselmodelabs($edgeselmodev)
    }
    $f.edgeselmode config -variable stldoctor.edgeselectmode
    $f.edgeselmode config -command { Ng_STLDoctor }
    global stldoctor.edgeselectmode
    pack $f.edgeselmode

    # edge buttons

    frame $f.edg -relief groove -borderwidth 3
    pack $f.edg -fill x

#    checkbutton $f.edg.bu -text "use external edges" \
#	-variable stldoctor.useexternaledges \
#	-command {Ng_STLDoctor; redraw}
#   pack $f.edg.bu -expand yes


    frame $f.edg.f0
    pack $f.edg.f0
    button $f.edg.f0.confirmedge -text "confirm" -command {Ng_STLDoctor confirmedge; redraw}
    button $f.edg.f0.candidateedge -text "candidate" -command {Ng_STLDoctor candidateedge; redraw}
    button $f.edg.f0.excludeedge -text "exclude" -command {Ng_STLDoctor excludeedge; redraw}
    button $f.edg.f0.undefinededge -text "undefined" -command {Ng_STLDoctor undefinededge; redraw}
    pack $f.edg.f0.confirmedge $f.edg.f0.candidateedge $f.edg.f0.excludeedge $f.edg.f0.undefinededge  -side left

    frame $f.edg.fa
    pack $f.edg.fa
    button $f.edg.fa.setallundefined -text "all undefined" -command {Ng_STLDoctor setallundefinededges; redraw}
    button $f.edg.fa.erasecandidates -text "candidates to undefined" -command {Ng_STLDoctor erasecandidateedges; redraw}
    pack $f.edg.fa.setallundefined $f.edg.fa.erasecandidates -side left


    frame $f.edg.fb
    pack $f.edg.fb
    button $f.edg.fb.confirmcandidates -text "candidates to confirmed" -command {Ng_STLDoctor confirmcandidateedges; redraw}
    button $f.edg.fb.confirmedtocandidates -text "confirmed to candidates" -command {Ng_STLDoctor confirmedtocandidateedges; redraw}
    pack $f.edg.fb.confirmcandidates $f.edg.fb.confirmedtocandidates -side left

    frame $f.edg.f1
    frame $f.edg.f2
    frame $f.edg.f3
    frame $f.edg.f4
    pack $f.edg.f1 $f.edg.f2 $f.edg.f3 $f.edg.f4

    button $f.edg.f1.exportedges -text "export edges" -command {Ng_STLDoctor exportedges}
    button $f.edg.f1.importedges -text "import edges" -command {Ng_STLDoctor importedges; redraw}
    button $f.edg.f1.saveedgedata -text "save edgedata" \
	-command { 
	    set types {
		{"Netgen Edgedata"   {.ned} } 
	    }
	    set file [tk_getSaveFile -filetypes $types -defaultextension ".ned"]
	    if {$file != ""} {
		Ng_STLDoctor saveedgedata $file
	}
    }

    button $f.edg.f1.loadedgedata -text "load edgedata" \
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

    button $f.edg.f1.importAVLedges -text "import AVL edges" \
	-command {
	    set types {{"Edge file"  {.edg }}}

	    set file [tk_getOpenFile -filetypes $types -defaultextension ".edg"]
	    if {$file != ""} {
		Ng_STLDoctor importexternaledges $file; 
	    }
	}

    pack $f.edg.f1.importAVLedges $f.edg.f1.loadedgedata $f.edg.f1.saveedgedata -side left

#    button $f.edg.f1.buildedges -text "build external edges" -command {Ng_STLDoctor buildexternaledges; redraw}
    frame $f.edg2 -relief groove -borderwidth 3
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
    set f [$wd.nb subwidget normals]

    frame $f.dt -relief groove -borderwidth 3
    pack $f.dt -fill x
    label $f.dt.lab -text "dirty triangle factor";
    entry $f.dt.ent -width 5 -relief sunken \
	-textvariable stldoctor.dirtytrigfact
    pack $f.dt.ent $f.dt.lab -side left  -expand yes

    frame $f.srt -relief groove -borderwidth 3
    pack $f.srt -fill x
    button $f.srt.bu -text "smooth reverted triangles geometric" -command {Ng_STLDoctor smoothrevertedtrigs; redraw }
    entry $f.srt.ent -width 5 -relief sunken \
	-textvariable stldoctor.smoothangle
    pack $f.srt.ent $f.srt.bu -side left  -expand yes

    frame $f.bdt -relief groove -borderwidth 3
    pack $f.bdt -fill x
    button $f.bdt.bu -text "mark dirty triangles" -command {Ng_STLDoctor markdirtytrigs; redraw }
    button $f.bdt.bu2 -text "smooth dirty triangles normal" -command {Ng_STLDoctor smoothdirtytrigs; redraw }
    pack $f.bdt.bu $f.bdt.bu2 -side left  -expand yes

    
    frame $f.sno -relief groove -borderwidth 3
    pack $f.sno
    
    label $f.sno.labrough -text "rough"
    scale $f.sno.scsmooth -orient horizontal -length 100 -from 0 -to 0.8 \
	-resolution 0.01 -variable stldoctor.smoothnormalsweight \
	-command { Ng_SetSTLParameters }
    label $f.sno.labsmooth -text "smooth"
    button $f.sno.smoothnormals -text "smooth normals" -command { Ng_STLDoctor smoothnormals; redraw}



    pack $f.sno.labrough $f.sno.scsmooth $f.sno.labsmooth $f.sno.smoothnormals -side left -padx 5

    frame $f.no -relief groove -borderwidth 3
    pack $f.no -fill x

    button $f.no.marknonsmoothnormals -text "mark non-smooth triangles" -command {Ng_STLDoctor marknonsmoothnormals; redraw}
    button $f.no.calcnormals -text "calculate normals from geometry" -command {Ng_STLDoctor calcnormals; redraw}

    pack  $f.no.marknonsmoothnormals $f.no.calcnormals -expand yes


    # ADVANCED **************************
    set f [$wd.nb subwidget advanced]


    frame $f.sc
    pack $f.sc -fill x
    checkbutton $f.sc.bu -text "spiral check" \
	-variable stldoctor.spiralcheck \
	-command {Ng_STLDoctor;}    
    checkbutton $f.sc.bu2 -text "cone check" \
	-variable stldoctor.conecheck \
	-command {Ng_STLDoctor;}    
    pack $f.sc.bu $f.sc.bu2


    tixControl $f.gtol -label "load-geometry tolerance factor" -integer false \
	-variable stldoctor.geom_tol_fact \
	-options {
	    entry.width 8
	    label.width 30
	    label.anchor e
	}	
    pack $f.gtol

    button $f.adap -text "Apply" -command {
	[.stldoctor_dlg.nb subwidget advanced].gtol invoke
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


	frame $w.vis -relief groove -borderwidth 3
	pack $w.vis

	checkbutton $w.vis.showfilledtrigs -text "Show filled triangles" \
	-variable viewoptions.drawfilledtrigs \
	-command { Ng_SetVisParameters; redraw }
	
	checkbutton $w.vis.showedges -text "Show edges" \
	-variable viewoptions.drawedges \
	-command { Ng_SetVisParameters; redraw }
	

	checkbutton $w.vis.showoutline -text "Show Triangle Outline" \
	-variable viewoptions.drawoutline \
	-command { Ng_SetVisParameters; redraw }

	pack $w.vis.showfilledtrigs  $w.vis.showoutline $w.vis.showedges

	tixControl $w.markedgedist -label "Mark edge dist: " -integer true \
	    -min 0 -max 999  \
	    -variable meshdoc.markedgedist \
	    -options {
		entry.width 3
		label.width 20
		label.anchor e
	    } \
	    -command {
		Ng_MeshDoctor markedgedist ${meshdoc.markedgedist}
		redraw
	    }
	pack $w.markedgedist
	
	button $w.deledge -text "Delete marked segments" -command {
	    Ng_MeshDoctor deletemarkedsegments
	    redraw
	}
	pack $w.deledge
	
	button $w.close -text "Close" -command { 
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



proc runtestdialog { } {
    source $::ngdir/ngshell.tcl
    set w .runtest_dlg
    
    if {[winfo exists .runtest_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w

	focus $w 
    } {
	toplevel $w

# in2d testing #
	frame $w.in2dframe 
	pack $w.in2dframe

        set in2dlogfile ""
 	tixLabelEntry $w.in2dframe.ent -label "in2d log-file: console if empty"  \
 	    -labelside top \
 	    -options {  
 		entry.textVariable in2dlogfile
 		entry.width 35
 		label.width 25
 		label.anchor w
 	    }	
 	button $w.in2dframe.btn -text "Browse" -command {
	    set types { { "Log file"   {.log}	} }
	    set in2dlogfile [tk_getOpenFile -filetypes $types -initialfile $in2dlogfile]
	}
 	button $w.in2dframe.test -text "Test in2d meshing" -command { ngtest in2d $in2dlogfile }

        
 	pack $w.in2dframe.test -side left -anchor s -padx 4 -pady 4
 	pack $w.in2dframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $w.in2dframe.btn -side left -anchor s -padx 4 -pady 4

        
# geo testing #
	frame $w.geoframe 
	pack $w.geoframe

        set geologfile "" 
 	tixLabelEntry $w.geoframe.ent -label "geo log-file: console if empty"  \
 	    -labelside top \
 	    -options {  
 		entry.textVariable geologfile
 		entry.width 35
 		label.width 25
 		label.anchor w
 	    }	
 	button $w.geoframe.btn -text "Browse" -command {
	    set types { { "Log file"   {.log}	} }
	    set geologfile [tk_getOpenFile -filetypes $types -initialfile $geologfile]
	}
 	button $w.geoframe.test -text "Test geo meshing" -command { ngtest geo $geologfile }

        
 	pack $w.geoframe.test -side left -anchor s -padx 4 -pady 4
 	pack $w.geoframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $w.geoframe.btn -side left -anchor s -padx 4 -pady 4

# stl testing #
	frame $w.stlframe 
	pack $w.stlframe

        set stllogfile ""
 	tixLabelEntry $w.stlframe.ent -label "stl log-file: console if empty"  \
 	    -labelside top \
 	    -options {  
 		entry.textVariable stllogfile
 		entry.width 35
 		label.width 25
 		label.anchor w
 	    }	
 	button $w.stlframe.btn -text "Browse" -command {
	    set types { { "Log file"   {.log}	} }
	    set stllogfile [tk_getOpenFile -filetypes $types -initialfile $stllogfile]
	}
 	button $w.stlframe.test -text "Test stl meshing" -command { ngtest stl $stllogfile }

        
 	pack $w.stlframe.test -side left -anchor s -padx 4 -pady 4
 	pack $w.stlframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $w.stlframe.btn -side left -anchor s -padx 4 -pady 4

# pde testing #
	frame $w.pdeframe 
	pack $w.pdeframe

        set pdelogfile ""
 	tixLabelEntry $w.pdeframe.ent -label "pde log-file: console if empty"  \
 	    -labelside top \
 	    -options {  
 		entry.textVariable pdelogfile
 		entry.width 35
 		label.width 25
 		label.anchor w
 	    }	
 	button $w.pdeframe.btn -text "Browse" -command {
	    set types { { "Log file"   {.log}	} }
	    set pdelogfile [tk_getOpenFile -filetypes $types -initialfile $pdelogfile]
	}
 	button $w.pdeframe.test -text "Test ngsolve pde's" -command { ngtest pde $pdelogfile }

        
 	pack $w.pdeframe.test -side left -anchor s -padx 4 -pady 4
 	pack $w.pdeframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4
 	pack $w.pdeframe.btn -side left -anchor s -padx 4 -pady 4
 
	wm title $w "Testing"
	focus .runtest_dlg 
    }
}

