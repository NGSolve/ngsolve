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
    #set hlist [$w.mtre subwidget hlist]

    set entities [Ng_GetOCCData getentities]
    set nrentities [expr [llength $entities]]


    if {$nrentities != 0} {
    
	#$hlist add Topology -itemtype text -text "Topology"
	
	#$hlist add Topology/CompSolids   -itemtype text -text "Composite Solids" -data "Composite Solids"
	#$hlist add Topology/FreeSolids   -itemtype text -text "Free Solids" -data "Free Solids"
	#$hlist add Topology/FreeShells   -itemtype text -text "Free Shells" -data "Free Shells"
	#$hlist add Topology/FreeFaces    -itemtype text -text "Free Faces" -data "Free Faces"
	#$hlist add Topology/FreeWires    -itemtype text -text "Free Wires" -data "Free Wires"
	#$hlist add Topology/FreeEdges    -itemtype text -text "Free Edges" -data "Free Edges"
	#$hlist add Topology/FreeVertices -itemtype text -text "Free Vertices" -data "Free Vertices"
    $w.tree insert {Topology} end -id "CompSolids" -text "Composite Solids"
    $w.tree insert {Topology} end -id "FreeSolids" -text "Free Solids"
    $w.tree insert {Topology} end -id "FreeShells" -text "Free Shells"
    $w.tree insert {Topology} end -id "FreeFaces" -text "Free Faces"
    $w.tree insert {Topology} end -id "FreeWires" -text "Free Wires"
    $w.tree insert {Topology} end -id "FreeEdges" -text "Free Edges"
    $w.tree insert {Topology} end -id "FreeVertices" -text "Free Vertices"
    
	#	$hlist add SingularEntities -itemtype text -text "Entities marked as singular"
    $w.tree item "Topology" -open true
	set i [expr 0]
	while {$i < $nrentities} {
	    set entity [lindex $entities [expr $i]]
	    incr i 1
	    set entityname [lindex $entities [expr $i]]
	    #$hlist add Topology/$entity -text $entityname -data $entityname
        set myroot [string range $entity 0 [string last / $entity]-1]
        $w.tree insert $myroot end -id $entity -text $entityname -value 1
	    incr i 1
	    #$w.mtre close Topology/$entity
	}
	
	#$w.mtre autosetmode
	
	#$w.mtre open Topology
	#$w.mtre close Topology/CompSolids
	#$w.mtre close Topology/FreeSolids
	#$w.mtre close Topology/FreeShells
	#$w.mtre close Topology/FreeFaces
	#$w.mtre close Topology/FreeWires
	#$w.mtre close Topology/FreeEdges
	#$w.mtre close Topology/FreeVertices
	
	set i [expr 0]
	while {$i < $nrentities} {
	    set entity [lindex $entities [expr $i]]
	    #$w.mtre close Topology/$entity
        $w.tree item $entity -open false
	    incr i 2
	}
	
	set faces [Ng_OCCCommand getunmeshedfaceinfo]    
	set nrfaces [expr [llength $faces]]
	if {$nrfaces >= 2} {
	    #$hlist add ErrorFaces -itemtype text -text "Faces with surface meshing error"
        $w.tree insert {} -id ErrorFaces -text "Faces with surface meshing error"
	    #$w.mtre open ErrorFaces
        $w.tree item ErrorFaces -open true
	    set i [expr 0]
	    while {$i < $nrfaces} {
		set entity [lindex $faces [expr $i]]
        set myroot [string range $entity 0 [string last / $entity]-1]
        if { [string length $myroot] == 0 } {
            set myroot ErrorFaces
        }
		incr i 1
		set entityname [lindex $faces [expr $i]]
		#$hlist add ErrorFaces/$entity -text $entityname -data $entityname
        $w.tree insert {myroot} end -id $entity -text $entityname -value 0
		incr i 1
	    }
	}
	

	set faces [Ng_OCCCommand getnotdrawablefaces]    
	set nrfaces [expr [llength $faces]]
	if {$nrfaces >= 2} {
	    #$hlist add NotDrawableFaces -itemtype text -text "Faces impossible to visualize"
        $w.tree insert {} -id NotDrawableFaces -text "Faces impossible to visualize"
	    #$w.mtre open NotDrawableFaces
        $w.tree item NotDrawableFaces -open true
	    set i [expr 0]
	    while {$i < $nrfaces} {
		set entity [lindex $faces [expr $i]]
        set myroot [string range $entity 0 [string last / $entity]-1]
        if { [string length $myroot ] == 0 } {
            set myroot NotDrawableFaces
        }
		incr i 1
		set entityname [lindex $faces [expr $i]]
		#$hlist add NotDrawableFaces/$entity -text $entityname -data $entityname
        $w.tree insert $myroot end -id $entity -text $entityname -value 0
		incr i 1
	    }
	}


	#$w.mtre autosetmode

	puts "done"
    }
}


proc rebuildoccdialog {} {
    if {[winfo exists .occ_dlg] == 1} {
    .occ_dlg.tree delete [.occ_dlg.tree children Topology]
	#[.occ_dlg.mtre subwidget hlist] delete all
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
	    # #set hlist [.occ_dlg.mtre subwidget hlist]
	    # #.occ_dlg.mtre open Topology
        # .occ_dlg.tree item Topology -open true
        # puts $entity2
	    # set slashpos [string last "/" $entity2]
	    # set entity3 [string range $entity2 0 [expr $slashpos-1]]
	    # while {$slashpos != -1} {
		# #.occ_dlg.mtre open Topology/$entity3
        # .occ_dlg.tree item $entity3 -open true
		# puts $entity3
		# set slashpos [string last "/" $entity3]
		# set entity3 [string range $entity3 0 [expr $slashpos-1]]
	    # }
	    #$hlist selection clear
	    .occ_dlg.tree see $entity2
	    .occ_dlg.tree selection set $entity2
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

        ttk::treeview $w.tree
        $w.tree insert {} end -id "Topology" -text "Topology"
        pack $w.tree -fill both -expand yes
	    occdialogbuildtree
	            
        bind $w.tree <Double-1> {
            set entityname [.occ_dlg.tree item [.occ_dlg.tree selection] -text ]
            set rootname ""
            if {[.occ_dlg.tree item [.occ_dlg.tree selection] -value] == 1 } {
                set rootname "Topology"
            }
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
	    ttk::button $w.cl -text "Close" -command {
		destroy .occ_dlg
	    }
	    
	    ttk::button $w.show -text "Show" -command {
        set entityname [.occ_dlg.tree item [.occ_dlg.tree selection] -text ]
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
	    ttk::button $w.hide -text "Hide" -command {
        set entityname [.occ_dlg.tree item [.occ_dlg.tree selection] -text ]        
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

	    ttk::button $w.swaporientation -text "Swap orientation" -command {
        set entityname [.occ_dlg.tree item [.occ_dlg.tree selection] -text ]                
		set spacepos [string first " " $entityname]
		set entitytype [string range $entityname 0 [expr $spacepos-1]]
		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		set spacepos2 [string first " " $helpstring]
		set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]

		Ng_OCCCommand swaporientation $entitytype $entitynumber
		set selectvisual geometry
		#	    Ng_SetVisParameters
		redraw
        .occ_dlg.tree delete [.occ_dlg.tree children Topology]
		#[.occ_dlg.mtre subwidget hlist] delete all
		occdialogbuildtree	
	    }

	    ttk::button $w.marksingular -text "Mark/Unmark as singular" -command {
        set entityname [.occ_dlg.tree item [.occ_dlg.tree selection] -text ]
		set spacepos [string first " " $entityname]
		if { $spacepos != 0 } {
		    set entitytype [string range $entityname 0 [expr $spacepos-1]]
		    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]
		    set spacepos2 [string first " " $helpstring]
		    if { $spacepos2 != 0 } {
			set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]
			puts $entitytype
            puts $entitynumber
			global ismarkedsingular
			Ng_OCCCommand marksingular $entitytype $entitynumber
						
			if { $ismarkedsingular == 0 } {
                .occ_dlg.tree tag remove "Color" [.occ_dlg.tree selection]
			} {
                .occ_dlg.tree tag add "Color" [.occ_dlg.tree selection]
                .occ_dlg.tree tag configure "Color" -foreground "red"
                .occ_dlg.tree tag configure "Color" -background "blue"
			}
		    }
		}                
	    }

	    ttk::checkbutton $w.zoomtohighlightedentity -text "Zoom to highlighted entity" \
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



	    ttk::frame $w.healing -relief groove -borderwidth 3

	    ttk::button $w.healing.checkentities -text "Analyze geometry" -command {
		set irregent [Ng_OCCCommand findsmallentities]
		set w .occ_dlg
		set nritems [expr [llength $irregent]]
		set i [expr 0]
        if {$nritems > 0 } {
            if { [.occ_dlg.tree exists ProblematicEntities] == 1 } {
                $w.tree delete ProblematicEntities
            }
            $w.tree insert {} end -id ProblematicEntities -text "Problematic Entities"
        }
		while {$i < $nritems} {
		    set entity [lindex $irregent [expr $i]]
		    incr i 1
		    set entityname [lindex $irregent [expr $i]]
            #puts $entity
            #puts $entityname
            set myroot [string range $entity 0 [string last / $entity]-1]
            if { [string length $myroot] == 0 } {
                set myroot ProblematicEntities
            }
            $w.tree insert $myroot end -id $entity -text $entityname            
		    incr i 1
		}
        $w.tree item ProblematicEntities -open true
	    }

	    # tixControl $w.healing.tolerance -label "Healing tolerance: " -integer false \
		# -variable occoptions.tolerance -min 1e-9 -max 1e6 \
		# -options {
		    # entry.width 6
		    # label.width 25
		    # label.anchor e
		# }	
        
        ttk::frame $w.healing.tolerance
        ttk::label $w.healing.tolerance.label -text "Healing tolerance: "
        ttk::spinbox $w.healing.tolerance.sp -textvariable occoptions.tolerance -width 6 -increment 0.01 -validate focus -validatecommand "my_validatespinbox %W %P 12" \
        -invalidcommand "my_invalidspinbox %W" -from -1e-9 -to 1e6 
        grid $w.healing.tolerance.label $w.healing.tolerance.sp
        
	    ttk::checkbutton $w.healing.fixsmalledges -text "Fix small edges" \
		-variable occoptions.fixsmalledges
	    
	    ttk::checkbutton $w.healing.fixspotstripfaces -text "Fix spot/strip faces" \
		-variable occoptions.fixspotstripfaces
	    
	    ttk::checkbutton $w.healing.sewfaces -text "Sew faces" \
		-variable occoptions.sewfaces
	    
	    ttk::checkbutton $w.healing.makesolids -text "Make solids" \
		-variable occoptions.makesolids
	    
	    ttk::checkbutton $w.healing.splitpartitions -text "Split partitions" \
		-variable occoptions.splitpartitions
	    
	    ttk::button $w.healing.heal -text "Heal geometry" -command { 
		#.occ_dlg.healing.tolerance.sp invoke
		Ng_OCCCommand shapehealing
		redraw 
        .occ_dlg.tree delete [.occ_dlg.tree children Topology]        
		#[.occ_dlg.mtre subwidget hlist] delete all        
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