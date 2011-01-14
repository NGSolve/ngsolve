.ngmenu.geometry add command -label "Scan CSG Geometry" -command { Ng_ParseGeometry }
.ngmenu.geometry add command -label "CSG Options..." -command geometryoptionsdialog

# only intern version !
# .ngmenu.geometry add separator
# .ngmenu.geometry add command -label "New Primitive"  \
#     -command newprimitivedialog -accelerator "<n><p>"
# .ngmenu.geometry add command -label "Edit Primitive" \
#     -command editprimitivedialog -accelerator "<e><p>"
# .ngmenu.geometry add command -label "Edit Solid" \
#     -command newsoliddialog -accelerator "<e><s>"
# .ngmenu.geometry add command -label "Choose Top Level " \
#     -command topleveldialog 
# .ngmenu.geometry add command -label "Identify" \
#     -command identifydialog 
.ngmenu.geometry add command -label "CSG Properties..." \
    -command topleveldialog2 











proc geometryoptionsdialog { } {


    set w .geometry_dlg
    
    if {[winfo exists .geometry_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {

	toplevel $w
	
	global geooptions
	
	Ng_GeometryOptions get

	checkbutton $w.drawcsg -text "Draw Geometry" \
	-variable geooptions.drawcsg 
	pack $w.drawcsg

	frame $w.fac
	pack $w.fac -pady 5
	label $w.fac.lab -text "Facets:";
	entry $w.fac.ent -width 8 -relief sunken \
	    -textvariable geooptions.facets
	pack $w.fac.lab $w.fac.ent  -side left
	
 
	frame $w.det
	pack $w.det -pady 5
	label $w.det.lab -text "Detail:";
	entry $w.det.ent -width 8 -relief sunken \
	    -textvariable geooptions.detail
	pack $w.det.lab $w.det.ent  -side left
	
	frame $w.cox
	pack $w.cox -pady 5
	label $w.cox.lab -text "min/max x:";
	entry $w.cox.ent1 -width 8 -relief sunken \
	    -textvariable geooptions.minx
	entry $w.cox.ent2 -width 8 -relief sunken \
	    -textvariable geooptions.maxx
	pack $w.cox.lab $w.cox.ent1 \
	    $w.cox.ent2  -side left
	
	frame $w.coy
	pack $w.coy -pady 5
	label $w.coy.lab -text "min/max y:";
	entry $w.coy.ent1 -width 8 -relief sunken \
	    -textvariable geooptions.miny
	entry $w.coy.ent2 -width 8 -relief sunken \
	    -textvariable geooptions.maxy
	pack $w.coy.lab $w.coy.ent1 \
	    $w.coy.ent2  -side left
	
	frame $w.coz
	pack $w.coz -pady 5
	label $w.coz.lab -text "min/max z:";
	entry $w.coz.ent1 -width 8 -relief sunken \
	    -textvariable geooptions.minz
	entry $w.coz.ent2 -width 8 -relief sunken \
	    -textvariable geooptions.maxz
	pack $w.coz.lab $w.coz.ent1 \
	    $w.coz.ent2  -side left
	


# 	tixButtonBox $w.bbox -orientation horizontal
# 	$w.bbox add ok    -text Apply  -underline 0 -width 5 \
# 	    -command { Ng_GeometryOptions set }
	
# 	$w.bbox add close -text Done   -underline 0 -width 5 \
# 	    -command { 
# 		Ng_GeometryOptions set
# 		destroy .geometry_dlg
# 	    }
# 	pack $w.bbox -side bottom -fill x
	

 	frame $w.bu
	pack $w.bu -fill x -ipady 3


 	button $w.bu.app -text "Apply" -command {
 	    Ng_GeometryOptions set
 	}
 	button $w.bu.ok -text "Done" -command {
 	    Ng_GeometryOptions set
 	    destroy .geometry_dlg
 	}
 	pack  $w.bu.app $w.bu.ok -side left -expand yes
    

	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Geometry options"
	focus $w
    }
}

















#
#
#   Edit primitive
#
#
proc editprimitivedialog2 { name } {

    global w classname

    set w .ep_dlg
    toplevel .$w

    Ng_GetPrimitiveData $name classname valuelist
        
    
    label $w.lab1 -text "Primitive Name:  $name";
    label $w.lab2 -text "Primitive Class: $classname";
    pack $w.lab1 $w.lab2 -fill x -pady 1m -padx 5m 
    
    frame $w.specific -relief groove

    global spec
    set spec(sphere) { cx cy cz rad }
    set spec(cylinder) { ax ay az bx by bz rad }
    set spec(plane) { px py pz nx ny nz }
    set spec(cone) { ax ay az bx by bz ra rb }
    set spec(brick) { p1x p1y p1z p2x p2y p2z p3x p3y p3z p4x p4y p4z } 
   
    set cnt 0
    foreach field $spec($classname) {

	frame $w.specific.f$cnt 
	pack $w.specific.f$cnt -side top -anchor ne

	label $w.specific.f$cnt.lab -text "$field"
	entry $w.specific.f$cnt.ent -textvariable dataval($cnt) \
	    -width 6 -relief sunken
	pack $w.specific.f$cnt.ent $w.specific.f$cnt.lab -side right
	$w.specific.f$cnt.ent delete 0 end
	$w.specific.f$cnt.ent insert 0 [lindex $valuelist $cnt]
	set cnt [expr $cnt + 1]
    }
    pack $w.specific


    button $w.cancel -text "cancel" -command {
	destroy $w 
    }

    button $w.ok -text "ok" -command {

	set valuelist ""
	set cnt 0
	foreach field $spec($classname) {
	    lappend valuelist $dataval($cnt)
	    set cnt [expr $cnt + 1]
	}
	Ng_SetPrimitiveData $name $valuelist
	destroy $w
    }
    pack  $w.cancel $w.ok -side left -expand yes

    bind $w <Return> { $w.ok  invoke}
    bind $w <Escape> { $w.cancel  invoke}
    

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w

#    grab $w
    focus $w.specific.f0.ent
}



#
#
#      Select primitve to edit
#
#

proc editprimitivedialog { } {
    global w

    set w .ep_dlg
    toplevel $w

    frame $w.frame -borderwidth 5m
    pack $w.frame -side top -expand yes -fill y

    listbox $w.frame.list -yscroll "$w.frame.scroll set" -setgrid 1 -height 12
    scrollbar $w.frame.scroll -command "$w.frame.list yview"
    pack $w.frame.scroll -side right -fill y
    pack $w.frame.list -side left -expand 1 -fill both
    

    Ng_GetPrimitiveList primlist
    foreach el $primlist {
	$w.frame.list insert end $el }

    button $w.cancel -text "cancel" -command { destroy $w }
    button $w.ok -text "ok" -command {
	set name [.ep_dlg.frame.list get active]
	puts "name=($name)"
	destroy $w
	if { $name != "" } { editprimitivedialog2 $name }
    }
    
    bind $w <Escape> { $w.cancel invoke }
    bind $w <Return> { $w.ok invoke }
    

    pack  $w.cancel $w.ok -side left -expand yes

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w

#    grab $w
    focus $w.frame.list
}



#
#
#         Create new primitive
#
#
proc newprimitivedialog { } {

    global w name

    set w .ap_dlg
    
    toplevel $w

    set name ""
    frame $w.f1
    pack $w.f1 -pady 2m
    label $w.f1.lab -text "Primitive Name: ";
    entry $w.f1.ent -width 5 -relief sunken \
	-textvariable name
    pack $w.f1.lab $w.f1.ent -side left
    
    frame $w.frame -borderwidth .5c
    pack $w.frame -side top -expand yes -fill y

    listbox $w.frame.list -yscroll "$w.frame.scroll set" -setgrid 1 -height 8 
    scrollbar $w.frame.scroll -command "$w.frame.list yview"
    pack $w.frame.scroll -side right -fill y
    pack $w.frame.list -side left -expand 1 -fill both
    
    $w.frame.list insert 0 sphere cylinder plane cone brick
    $w.frame.list activate 0
    
    button $w.ok -text "ok" -command {
	Ng_CreatePrimitive [$w.frame.list get active]  $name
	destroy $w
	editprimitivedialog2 $name
    }

    button $w.cancel -text "cancel" -command {
	destroy $w
    }
    
    pack  $w.cancel $w.ok -side left -expand yes -pady 2m


    bind $w <Escape> { $w.cancel invoke }
    bind $w <Return> { $w.ok invoke }

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w

#    grab $w
    focus $w.f1.ent
}






proc newsoliddialog { } {

    global w name val sollist

    set w .ns_dlg
    toplevel $w

    set name ""
    frame $w.f1
    label $w.f1.lab -text "Solid Name: ";
    entry $w.f1.ent -width 5 -relief sunken \
	-textvariable name
    $w.f1.ent delete 0 end
    button $w.f1.getsel -text "Get Selected" -command { 
	$w.f1.ent delete 0 end
	$w.f1.ent insert 0 [$w.f3.list get active]
	$w.bu.get invoke
    }
    pack $w.f1.getsel -side bottom
    pack $w.f1.ent $w.f1.lab -side right


    frame $w.f3 -borderwidth .5c
    listbox $w.f3.list -yscroll "$w.f3.scroll set" -setgrid 1 -height 12
    scrollbar $w.f3.scroll -command "$w.f3.list yview"
    pack $w.f3.scroll -side right -fill y
    pack $w.f3.list -side left -expand 1 -fill both
    
    Ng_GetSolidList sollist
    foreach el $sollist {
	$w.f3.list insert end $el }

    frame $w.f2
    label $w.f2.lab -text "Solid Description: ";
    pack $w.f2.lab


    entry $w.f2.ent -width 100 -relief sunken \
	-textvariable val  -xscrollcommand "$w.f2.scr set"
    scrollbar $w.f2.scr -relief sunken -orient horiz -command \
	"$w.f2.ent xview"
    $w.f2.ent delete 0 end
    pack $w.f2.ent $w.f2.scr -fill x



    frame $w.bu
    button $w.bu.close -text "close" -command {
	destroy $w
    }

    button $w.bu.get -text "get data" -command {
	Ng_GetSolidData $name val
    }

    button $w.bu.set -text "set data" -command {
	Ng_SetSolidData $name $val
    }

    pack $w.bu.get $w.bu.set $w.bu.close -side left 


    pack $w.bu -pady 5 -side bottom                     ;# buttons
    pack $w.f2 -pady 5 -side bottom                     ;# edit field
    pack $w.f1 -pady 5 -side left                       ;# name
    pack $w.f3 -side left -expand yes -fill y           ;# listbox



    bind $w <Escape> { $w.bu.close invoke }

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w

#    grab $w
    focus $w
}







#
#  Edit top level objects
#
#


proc toplevelproperties { w solname surfname } {

    global properties

    Ng_TopLevel getprop $solname $surfname properties


    set w .tlprop_dlg

    if {[winfo exists $w] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w
    
	label $w.lab1 -text "Red"
	scale $w.scale1 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(red)
	
	label $w.lab2 -text "Green"
	scale $w.scale2 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	-command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(green)
	
	label $w.lab3 -text "Blue"
	scale $w.scale3 -orient horizontal -length 300 -from 0 -to 1 \
	    -resolution 0.01  -tickinterval 0.2 \
	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(blue)

	
	pack $w.lab1 $w.scale1 $w.lab2 $w.scale2 $w.lab3 $w.scale3

	checkbutton $w.cb4 -text "Visible" \
	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } \
	    -variable properties(visible)
	
	checkbutton $w.cb5 -text "Transparent" \
	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } \
	    -variable properties(transp)
	
	
	pack  $w.cb4 $w.cb5
	
	
	frame $w.bu
	pack $w.bu -fill x
	button $w.bu.ok -text "Ok" -command "destroy .tlprop_dlg"
	pack $w.bu.ok  -expand yes
    
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	focus $w
    }
    wm title $w "Properties $solname $surfname"
}






proc topleveldialog { } {

    global w name val sollist

    set w .tl_dlg
    toplevel $w



    frame $w.sol -borderwidth .5c
    listbox $w.sol.list -yscroll "$w.sol.scroll set" -setgrid 1 -height 12
    scrollbar $w.sol.scroll -command "$w.sol.list yview"
    pack $w.sol.scroll -side right -fill y
    pack $w.sol.list -side left -expand 1 -fill both
    
    Ng_GetSolidList sollist
    foreach el $sollist {
	$w.sol.list insert end $el }
    Ng_GetPrimitiveList sollist
    foreach el $sollist {
	$w.sol.list insert end $el }




    frame $w.sul -borderwidth .5c
    listbox $w.sul.list -yscroll "$w.sul.scroll set" -setgrid 1 -height 12
    scrollbar $w.sul.scroll -command "$w.sul.list yview"
    pack $w.sul.scroll -side right -fill y
    pack $w.sul.list -side left -expand 1 -fill both
    
    Ng_GetSurfaceList sollist
    foreach el $sollist {
	$w.sul.list insert end $el }





    frame $w.topl -borderwidth .5c
    listbox $w.topl.list -yscroll "$w.topl.scroll set" -setgrid 1 -height 12 \
	-command { puts hi }
    scrollbar $w.topl.scroll -command "$w.topl.list yview"
    pack $w.topl.scroll -side right -fill y
    pack $w.topl.list -side left -expand 1 -fill both
    
    Ng_TopLevel getlist sollist
    puts $sollist
    foreach el $sollist {
	set hel "[ lindex $el 0 ]"
	if { [ llength $el ] == 2 } {
	    set hel "[ lindex $el 1 ] on [ lindex $el 0 ]"
	}
	$w.topl.list insert end $hel 
    }


    frame $w.bu

    button $w.bu.close -text "close" -command {
	destroy $w
    }
    button $w.bu.addsol -text "Add Solid" -command {
	set solname [$w.sol.list get active]
	Ng_TopLevel set $solname ""
	Ng_ParseGeometry
	$w.topl.list insert end $solname
    }

    button $w.bu.addsurf -text "Add Surface" -command {
	set solname [$w.sol.list get active]
	set surfname [$w.sul.list get active]
	Ng_TopLevel set $solname $surfname
	Ng_ParseGeometry
	puts "$solname on $surfname"
	$w.topl.list insert end "$surfname on $solname"
    }

    button $w.bu.remsol -text "Remove" -command {
	set solname [$w.topl.list get active]
	set surfname ""
	if { [llength $solname] == 3 } {
	    set surfname [lindex $solname 0]
	    set solname [lindex $solname 2]
	}
	Ng_TopLevel remove $solname $surfname
	Ng_ParseGeometry
	$w.topl.list delete active
    }

    button $w.bu.prop -text "Properties" -command {
	set solname [$w.topl.list get active]
	set surfname ""
	if { [llength $solname] == 3 } {
	    set surfname [lindex $solname 0]
	    set solname [lindex $solname 2]
	}
	toplevelproperties tlp $solname $surfname
    }


    

    pack  $w.bu.close $w.bu.addsol $w.bu.addsurf $w.bu.remsol $w.bu.prop  -side left 


    pack $w.bu -side bottom
    pack $w.sol -side left -expand yes -fill y           ;# listbox
    pack $w.sul -side left -expand yes -fill y           ;# listbox
    pack $w.topl -side left -expand yes -fill y           ;# listbox


    bind $w <Escape> { $w.bu.close invoke }

    wm withdraw $w
    wm geom $w +100+100
    wm deiconify $w

#    grab $w
    focus $w
}






proc topleveldialog2 { } {
    set w .tl2_dlg
    
    if {[winfo exists .tl2_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

	global name val sollist

	frame $w.topl -borderwidth .5c
	listbox $w.topl.list -yscroll "$w.topl.scroll set" -setgrid 1 -height 12
	scrollbar $w.topl.scroll -command "$w.topl.list yview"
	pack $w.topl.scroll -side right -fill y
	pack $w.topl.list -side left -expand 1 -fill both
	
	Ng_TopLevel getlist sollist
	puts $sollist
	set i 1
	foreach el $sollist {
	    set hel "$i: [ lindex $el 0 ]"
	    if { [ llength $el ] == 2 } {
		set hel "$i: [ lindex $el 1 ] on [ lindex $el 0 ]"
	    }
	    incr i
	    $w.topl.list insert end $hel }
	
	
	frame $w.bu
	
	button $w.bu.close -text "close" -command {
	    destroy .tl2_dlg
	}
	

	button $w.bu.prop -text "Properties" -command {
	    set solname [.tl2_dlg.topl.list get active]
	    set surfname ""
	    if { [llength $solname] == 2 } {
		set solname [lindex $solname 1]
	    }
	    if { [llength $solname] == 4 } {
		set surfname [lindex $solname 1]
		set solname [lindex $solname 3]
	    }
	    toplevelproperties tlp $solname $surfname
	}
	
	pack $w.bu.close $w.bu.prop  -side left 
	pack $w.bu -side bottom
	pack $w.topl -side left -expand yes -fill y           ;# listbox
		
	bind .tl2_dlg.topl.list <Double-1> {
	    set solname [.tl2_dlg.topl.list get  @%x,%y]
	    set surfname ""
	    if { [llength $solname] == 2 } {
		set solname [lindex $solname 1]
	    }
	    if { [llength $solname] == 4 } {
		set surfname [lindex $solname 1]
		set solname [lindex $solname 3]
	    }
	    toplevelproperties tlp $solname $surfname
	}
	
	bind .tl2_dlg <Escape> { .tl2_dlg.bu.close invoke }
	
	wm withdraw $w
	wm geom $w +100+100
	wm deiconify $w
	wm title $w "Top-Level Options"	
	focus $w
    }
}





