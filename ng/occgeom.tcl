if { [catch { load liboccvis[info sharedlibextension] Ng_OCC } result ] } {
    # puts "cannot load occ" 
    # puts "error: $result"
}


if { [catch { Ng_OCCCommand isoccgeometryloaded }] } {
    # dummy 	
    proc rebuildoccdialog { } { }
} {
    puts "OCC module loaded"
    set hasocc yes



.ngmenu.geometry add separator

.ngmenu.geometry add command -label "IGES/STEP Topology Explorer/Doctor..." \
    -command { occdialog; }


# Philippose - 30/01/2009
# Add menu item for local face mesh size definition in the 
# TCL Gui
.ngmenu.geometry add command -label "Edit Face Mesh Size..." \
    -command { surfacemeshsizedialog }


# .ngmenu.geometry add command -label "OCC Construction" \
    #     -command { Ng_OCCConstruction; }




set entities [ ]


proc occdialogbuildtree {} {
    global entities

    set w .occ_dlg
    set hlist [$w.mtre subwidget hlist]

    set entities [Ng_GetOCCData getentities]
    set nrentities [expr [llength $entities]]


    if {$nrentities != 0} {

	$hlist add Topology -itemtype text -text "Topology"
	
	$hlist add Topology/CompSolids   -itemtype text -text "Composite Solids" -data "Composite Solids"
	$hlist add Topology/FreeSolids   -itemtype text -text "Free Solids" -data "Free Solids"
	$hlist add Topology/FreeShells   -itemtype text -text "Free Shells" -data "Free Shells"
	$hlist add Topology/FreeFaces    -itemtype text -text "Free Faces" -data "Free Faces"
	$hlist add Topology/FreeWires    -itemtype text -text "Free Wires" -data "Free Wires"
	$hlist add Topology/FreeEdges    -itemtype text -text "Free Edges" -data "Free Edges"
	$hlist add Topology/FreeVertices -itemtype text -text "Free Vertices" -data "Free Vertices"

	#	$hlist add SingularEntities -itemtype text -text "Entities marked as singular"
	
	set i [expr 0]
	while {$i < $nrentities} {
	    set entity [lindex $entities [expr $i]]
	    incr i 1
	    set entityname [lindex $entities [expr $i]]
	    $hlist add Topology/$entity -text $entityname -data $entityname
	    incr i 1
	    $w.mtre close Topology/$entity
	}
	
	$w.mtre autosetmode
	
	$w.mtre open Topology
	$w.mtre close Topology/CompSolids
	$w.mtre close Topology/FreeSolids
	$w.mtre close Topology/FreeShells
	$w.mtre close Topology/FreeFaces
	$w.mtre close Topology/FreeWires
	$w.mtre close Topology/FreeEdges
	$w.mtre close Topology/FreeVertices
	
	set i [expr 0]
	while {$i < $nrentities} {
	    set entity [lindex $entities [expr $i]]
	    $w.mtre close Topology/$entity
	    incr i 2
	}
	
	set faces [Ng_OCCCommand getunmeshedfaceinfo]    
	set nrfaces [expr [llength $faces]]
	if {$nrfaces >= 2} {
	    $hlist add ErrorFaces -itemtype text -text "Faces with surface meshing error"
	    $w.mtre open ErrorFaces
	    set i [expr 0]
	    while {$i < $nrfaces} {
		set entity [lindex $faces [expr $i]]
		incr i 1
		set entityname [lindex $faces [expr $i]]
		$hlist add ErrorFaces/$entity -text $entityname -data $entityname
		incr i 1
	    }
	}
	

	set faces [Ng_OCCCommand getnotdrawablefaces]    
	set nrfaces [expr [llength $faces]]
	if {$nrfaces >= 2} {
	    $hlist add NotDrawableFaces -itemtype text -text "Faces impossible to visualize"
	    $w.mtre open NotDrawableFaces
	    set i [expr 0]
	    while {$i < $nrfaces} {
		set entity [lindex $faces [expr $i]]
		incr i 1
		set entityname [lindex $faces [expr $i]]
		$hlist add NotDrawableFaces/$entity -text $entityname -data $entityname
		incr i 1
	    }
	}


	$w.mtre autosetmode

	puts "done"
    }
}


proc rebuildoccdialog {} {
    if {[winfo exists .occ_dlg] == 1} {
	[.occ_dlg.mtre subwidget hlist] delete all
	occdialogbuildtree 
    }
}

proc checkoccloaded { } {
    set isoccgeometryloaded [Ng_OCCCommand isoccgeometryloaded]
    if {$isoccgeometryloaded == 0} {
	puts "no IGES/STEP geometry loaded"
	destroy .occ_dlg
    }
}

#proc setocctolerance { } {
#    set w .setocctolerance
#}


proc selectentity { entityname } {
    global entities
    set nrentities [expr [llength $entities]]
    set i [expr 0]
    while {$i < $nrentities} {
	set entitylength []
	
	set entity2 [lindex $entities [expr $i]]
	incr i 1
	set entityname2 [lindex $entities [expr $i]]
	incr i 1
	set entityname2 [string range $entityname2 0 [expr [string length $entityname]-1]]
	
	if {$entityname == $entityname2} {
	    set hlist [.occ_dlg.mtre subwidget hlist]
	    .occ_dlg.mtre open Topology
	    set slashpos [string last "/" $entity2]
	    set entity3 [string range $entity2 0 [expr $slashpos-1]]
	    while {$slashpos != -1} {
		.occ_dlg.mtre open Topology/$entity3
		
		set slashpos [string last "/" $entity3]
		set entity3 [string range $entity3 0 [expr $slashpos-1]]
	    }
	    $hlist selection clear
	    $hlist see Topology/$entity2
	    $hlist selection set Topology/$entity2
	} 
    }	    
}



proc occdialog { } {
    
    uplevel 1 {

	global entities
	set selectvisual geometry
	Ng_SetVisParameters
	redraw

	set w .occ_dlg
	
	if {[winfo exists .occ_dlg] == 1} {
	    wm withdraw $w
	    wm deiconify $w
	    focus $w 
	} {	
	    toplevel $w

	    tixTree $w.mtre -options { separator "/" }
	    pack $w.mtre -fill both -expand yes

	    occdialogbuildtree

	    set hlist [$w.mtre subwidget hlist]
            $hlist configure -selectforeground black
            $hlist configure -selectbackground grey

	    set solname {""}

	    
	    bind $hlist <Double-1> {
		set oldsolname {$solname}
		set solname [[.occ_dlg.mtre subwidget hlist] info selection]
		if {$solname != "" && $oldsolname != $solname } {
		    set seppos [string first "/" $solname]
		    set rootname [string range $solname 0 [expr $seppos-1]]

		    set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]
		    set spacepos [string first " " $entityname]
		    set entitytype [string range $entityname 0 [expr $spacepos-1]]
		    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		    set spacepos2 [string first " " $helpstring]
		    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]
		    if {$rootname == "Topology"} {
			Ng_OCCCommand highlightentity $entitytype $entitynumber
			set selectvisual geometry
			redraw
		    } {
			set brackpos [string first " (" $entityname]
			if {$brackpos != -1} {
			    set entityname [string range $entityname 0 $brackpos]
			}

			selectentity $entityname
		    }
		}
	    }
	    
	    button $w.cl -text "Close" -command {
		destroy .occ_dlg
	    }
	    
	    button $w.show -text "Show" -command {
		set solname [[.occ_dlg.mtre subwidget hlist] info selection]
		set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]
		set spacepos [string first " " $entityname]
		set entitytype [string range $entityname 0 [expr $spacepos-1]]
		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		set spacepos2 [string first " " $helpstring]
		set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]

		Ng_OCCCommand show $entitytype $entitynumber
		set selectvisual geometry
		#	    Ng_SetVisParameters
		redraw
	    }
	    button $w.hide -text "Hide" -command {
		set solname [[.occ_dlg.mtre subwidget hlist] info selection]
		set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]
		set spacepos [string first " " $entityname]
		set entitytype [string range $entityname 0 [expr $spacepos-1]]
		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		set spacepos2 [string first " " $helpstring]
		set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]

		Ng_OCCCommand hide $entitytype $entitynumber
		set selectvisual geometry
		#	    Ng_SetVisParameters
		redraw
	    }

	    button $w.swaporientation -text "Swap orientation" -command {
		set solname [[.occ_dlg.mtre subwidget hlist] info selection]
		set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]
		set spacepos [string first " " $entityname]
		set entitytype [string range $entityname 0 [expr $spacepos-1]]
		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		set spacepos2 [string first " " $helpstring]
		set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]

		Ng_OCCCommand swaporientation $entitytype $entitynumber
		set selectvisual geometry
		#	    Ng_SetVisParameters
		redraw

		[.occ_dlg.mtre subwidget hlist] delete all
		occdialogbuildtree	
	    }

	    button $w.marksingular -text "Mark/Unmark as singular" -command {
		set solname [[.occ_dlg.mtre subwidget hlist] info selection]
		set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]
		set spacepos [string first " " $entityname]
		if { $spacepos != 0 } {
		    set entitytype [string range $entityname 0 [expr $spacepos-1]]
		    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		    set spacepos2 [string first " " $helpstring]
		    if { $spacepos2 != 0 } {
			set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]
			
			global ismarkedsingular
			Ng_OCCCommand marksingular $entitytype $entitynumber
			
			set hlist [$w.mtre subwidget hlist]
			
			#	    $hlist entryconfigure $solname -text "hallo"
			#	    set style1 [tixDisplayStyle imagetext -font 8x13]
			set style1 [tixDisplayStyle imagetext -foreground black -background white -selectforeground white -selectbackground blue]
			set style2 [tixDisplayStyle imagetext -foreground red -background white -selectforeground red -selectbackground blue]
			
			if { $ismarkedsingular == 0 } {
			    $hlist entryconfigure $solname -style $style1
			} {
			    $hlist entryconfigure $solname -style $style2
			}

			#		    set hlist [$w.mtre subwidget hlist]
			#		    foreach solname2 $hlist {
			#			if { $ismarkedsingular == 0 } {
			#			    $hlist entryconfigure $solname2 -style $style1
			#			} {
			#			    $hlist entryconfigure $solname2 -style $style2
			#			}
			#		    }
		    }
		}
		#	    $hlist add test -after $solname

		#	    $hlist add SingularEntities/$entityname -text $entityname
		#	    set selectvisual geometry
		#	    Ng_SetVisParameters
		#	    redraw

		#	    [.occ_dlg.mtre subwidget hlist] delete all
		#	    occdialogbuildtree	
	    }


	    checkbutton $w.zoomtohighlightedentity -text "Zoom to highlighted entity" \
		-variable occoptions.zoomtohighlightedentity \
		-command {
		    Ng_SetOCCVisParameters
		    if { ${occoptions.zoomtohighlightedentity} == 1} {
			set selectvisual geometry
			#		    Ng_SetVisParameters
			Ng_OCCCommand redrawstatus 1
			redraw
		    } {
			Ng_OCCCommand redrawstatus 0
		    }
		}



	    frame $w.healing -relief groove -borderwidth 3

	    button $w.healing.checkentities -text "Analyze geometry" -command {
		set irregent [Ng_OCCCommand findsmallentities]

		set w .occ_dlg
		set hlist [$w.mtre subwidget hlist]
		
		$hlist add ProblematicEntities -text "Problematic Entities"
		$hlist delete offsprings ProblematicEntities

		set nritems [expr [llength $irregent]]
		set i [expr 0]
		while {$i < $nritems} {
		    set entity [lindex $irregent [expr $i]]
		    incr i 1
		    set entityname [lindex $irregent [expr $i]]
		    $hlist add ProblematicEntities/$entity -text $entityname -data $entityname
		    incr i 1
		}
		$w.mtre open ProblematicEntities
		$w.mtre autosetmode
	    }

	    tixControl $w.healing.tolerance -label "Healing tolerance: " -integer false \
		-variable occoptions.tolerance -min 1e-9 -max 1e6 \
		-options {
		    entry.width 6
		    label.width 25
		    label.anchor e
		}	

	    checkbutton $w.healing.fixsmalledges -text "Fix small edges" \
		-variable occoptions.fixsmalledges
	    
	    checkbutton $w.healing.fixspotstripfaces -text "Fix spot/strip faces" \
		-variable occoptions.fixspotstripfaces
	    
	    checkbutton $w.healing.sewfaces -text "Sew faces" \
		-variable occoptions.sewfaces
	    
	    checkbutton $w.healing.makesolids -text "Make solids" \
		-variable occoptions.makesolids
	    
	    checkbutton $w.healing.splitpartitions -text "Split partitions" \
		-variable occoptions.splitpartitions
	    
	    button $w.healing.heal -text "Heal geometry" -command { 
		.occ_dlg.healing.tolerance invoke
		Ng_OCCCommand shapehealing
		redraw 
		[.occ_dlg.mtre subwidget hlist] delete all
		occdialogbuildtree
	    }

	    pack $w.healing.checkentities

	    pack $w.healing.tolerance $w.healing.fixsmalledges \
		$w.healing.fixspotstripfaces $w.healing.sewfaces \
		$w.healing.makesolids $w.healing.splitpartitions -anchor w

	    pack $w.healing.heal	




	    pack $w.show $w.hide

	    pack $w.zoomtohighlightedentity -anchor w
	    #	pack $w.checkentities
	    pack $w.swaporientation
	    pack $w.marksingular
	    pack $w.healing -fill x
	    pack $w.cl
	    
	    
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "IGES/STEP Topology Explorer/Doctor"
	    focus .occ_dlg
	}
    }
}


}