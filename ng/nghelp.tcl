proc print_commandline_help { } {
    
    puts "Usage: ng { options }"

    puts "-geofile=filename      Input geometry file (alternative:  ng filename)"
    puts "-meshfile=filename     Output mesh file"
    puts "-verycoarse, -coarse, -moderate, -fine, -veryfine"
    puts "                       Automatic mesh-size selection"
    puts "-meshsizefile=filename Load mesh-size file with local mesh sizes"
    puts "-meshfiletype={\"Neutral Format\", ...}"
    puts "                       Filetype of output file, default is netgen file"
    puts "-batchmode             Run Netgen in batchmode"
    puts "-inputmeshfile=filename"
    puts "                       Input mesh file (batchmode only)"
    puts "-mergefile=filename    Merge with mesh file (batchmode only)"
    puts "-refinementfile=filename"
    puts "                       Use refinementinfo from file (batchmode only)"
    puts "-serversocket=\#num    Start a Netgen server with port \#num"
    puts "-V                     Print additional information"
    puts "-testout=filename      file for test output"

    if { [catch { NGS_GetData } ] == 0 } { 
	puts "\nNGSolve parameters:"
	puts "-pdefile=filename      Load pde input file"
	puts "-solve                 Solve pde once"
	puts "-solve=n               Solve pde by n adaptive refinement steps"
	puts "-recent                Load and solve most recently loaded pde"
    }

}



proc set_menu_help { entry helpmsg } {
    global menuhelps
    set menuhelps($entry) $helpmsg
}

proc show_menu_help { entry } {
    global menuhelps


    if {[catch {set helptext $menuhelps($entry)}]} {
	set helptext "no help available   "
    }    

    .helpline configure -text $helptext
    
    if {[winfo exists .senshelp_dlg]==1} {
	.senshelp_dlg.text delete 1.0 end
	.senshelp_dlg.text insert end "Menu item: $entry\n\n"
	.senshelp_dlg.text insert end $helptext
    }
}


#tixBalloon .balloon -statusbar .helpline

proc set_control_help { control helpmsg } {
    bind $control <Enter> "show_control_help {$helpmsg}"
    bind $control <Leave> "show_control_help {None}"
    .balloon bind  $control -balloonmsg $helpmsg -statusmsg $helpmsg
#    puts "Add Help to $control"
}

proc show_control_help { helpmsg } {
    .helpline configure -text $helpmsg
    if {[winfo exists .senshelp_dlg]==1} {
	.senshelp_dlg.text delete 1.0 end
	.senshelp_dlg.text insert end $helpmsg
    }
}


proc sensitivehelpdialog { show } {

    set w .senshelp_dlg
    
    if {[winfo exists .senshelp_dlg] == 1} {

	if { $show == 1 } {
	    wm withdraw .senshelp_dlg
	    wm deiconify $w
	    focus $w 
	} {
	    wm withdraw $w
	}
    } {
	toplevel $w
#	wm minsize $w 200 150

	global senshelptext

	text $w.text -yscrollcommand "$w.scroll set" -setgrid true \
	    -width 40 -height 10  -wrap word
	scrollbar $w.scroll -command "$w.text yview"
	pack $w.scroll -side right -fill y
	pack $w.text -expand yes -fill both

	frame $w.bu
	pack $w.bu
	# -fill x
	
	button $w.close -text "Close" \
	    -command { 
		wm withdraw .senshelp_dlg
		set showsensitivehelp 0
	    }
	pack $w.close
    
	
	if { $show == 1 } {
	    wm withdraw $w
	    wm geom $w +100+100
	    wm deiconify $w
	    wm title $w "Help"
	    focus $w
	}
    }
}



set_menu_help "File"  "In File menu you can load and store geometries, meshes etc." 

set_menu_help "New Geometry"  "Deletes current geometry"
set_menu_help "Load Geometry"  "Loads Geometry file in one of the formats STL (ASCII or binary), Constructive Solid Geometry (.geo) or 2D geometry. Please have a look into Netgen User's manuel for more details."
set_menu_help "Save Geometry" "Saves STL Geometry in in either ASCII or binary STL format."
set_menu_help "Load Mesh" "Loads surface and volume mesh in Netgen internal format."
set_menu_help "Save Mesh" "Saves surface and volume mesh in Netgen internal format."
set_menu_help "Write EPS File" "Dumps OpenGL rendering to EPS File."
set_menu_help "Save Options" "Saves current options in file \"ng.opt\". These options will be loaded again when starting ng in the same directory."
set_menu_help "Export Mesh" "Exports mesh in format defined by Export Filetype."
set_menu_help "Export Filetype" "Selects file format for exporting mesh. Please have a look into the Netgen User's manual for more information."
set_menu_help "Import Mesh" "Imports surface or volume mesh in exchange format."
set_menu_help "Quit" "Quits Netgen"

set_menu_help "Geometry" "Preparing geometries, visualiztion of geometries."
set_menu_help "Scan CSG Geometry" "Generates surface triangulation for rendering"
set_menu_help "CSG Options" "Sets Options for CSG visualization (bounding box, detail size, number of facets)."
set_menu_help "CSG Properties" "Defines appearence of current CSG geometry (color, visibility, transparency)"
set_menu_help "STL Doctor" "Calls STL Doctor for preprocessing STL geometry files."
set_menu_help "STL Info" "Retrieves information about current STL geometry."

set_menu_help "Mesh" "Menu for mesh generation"
set_menu_help "Generate Mesh" "Generates mesh from geometry, same as Button \"Generate Mesh\""
set_menu_help "Stop Meshing" "Terminates meshgeneration. It may take a while until meshing terminates, please be patient."
set_menu_help "Meshing Options" "Set options for mesh generation."
set_menu_help "Delete Mesh" "Deletes mesh. Not necessary before generation of new mesh."
set_menu_help "Delete Vol Mesh" "Deletes only volume mesh."
set_menu_help "Mesh Quality" "Computs element shape measures. Triangle angles are inner angles of all triangles (faces of tetrahedra). Tet angles are angles between faces of tetrahedra."
set_menu_help "Check Surface Mesh" "Checks consistency and overlap of surface mesh. Marks overlapping elements as bad elements, please enable visualization of bad elements in View->Mesh."
set_menu_help "Check Volume Mesh" "Checks conformity of volume mesh."
set_menu_help "Edit Boundary Conditions" "Open dialog for setting boundary condition numbers for individual faces."
set_menu_help "Analyze Geometry" "Perform only first step in mesh generation. Action depends on geometry type, e.g. generates charts for STL mesh, find vertices in CSG geometries."
set_menu_help "Mesh Edges" "Meshes edges"
set_menu_help "Mesh Surface" "Generates surface mesh. Includes already surface optimization for some geomtry types."
set_menu_help "Optimize Surface" "Optimizes surface mesh."
set_menu_help "Surface Optim. Step" "Performs a specific surface optimiztion step. Mesh smoothing moves nodes. edge swapping swaps the diagonal of a quadrilateral built by two triangles, criterion either by number of nodes, or anlges. Combine points eliminates triangles by combining points (in the center of gravity)."
set_menu_help "Mesh Volume" "Performs volume meshing. Algorithm is a combination of Delaunay and Rule-based Advancing Front"
set_menu_help "Optimize Volume" "Performs additional volume optimization steps"
set_menu_help "Smooth Opt Volume" "Performs optimization steps by smoothing iterations"
set_menu_help "Smooth Opt Volume Jacobian" "Volume optimization by smoothing iterations. Criterion is optimization of Jacobi determinants. This optimization step is also available for 10-node tetrahedra."

set_menu_help "View" "Sets viewing options"
set_menu_help "Zoom all" "Zooms scene to show whole object"
set_menu_help "Center" "Defines center of rotation"
set_menu_help "Viewing Options" "Sets viewing options for geometry, mesh, lighting"
set_menu_help "Clipping Plane" "Introduces clipping plane. The clipping plane is defined by the normal vector, and a scaled offset. Clipping of performed by OpenGl rendering"
set_menu_help "Quality Plot" "Shows the element quality distribution histogram. Measure is volume scaled by edge-length to the third. Optimal elements have measure 1."
set_menu_help "Sensitve Help" "Shows this help window"

set_menu_help "Mesh-size" "Manipulations of existing mesh"
set_menu_help "Refine uniform" "Refines mesh by splitting elements into eight childs (algorithm of J. Bey)"
set_menu_help "Second Order" "Converts 4 node elements to 10 node elements. Edge-midpoitns are projected to the geometry."
set_menu_help "Refinement Dialog" "Controls local mesh refinement"
set_menu_help "Load Meshsize" "Loads mesh-size file for local mesh refinement."
set_menu_help "MS from Surf Mesh" "Defines mesh-size by the surface mesh."




set f .options_dlg.nb.nbframe.general
# set_control_help $f "General meshing page"
set_control_help $f.fine "Controls relative mesh size.\nThis control affects other mesh-size controls in common"
set_control_help $f.first "First step in mesh generation. Usually, meshing  starts from \"analyze geometry\". If the surface mesh is already available \"First step\" should be set to \"mesh volume\""
set_control_help $f.last "Last step in mesh generation. If only the surface mesh is required, please set \"Last Step\" to \"Optimize Surface\""

set_control_help .bubar.surfm "Start mesh generation"
set_control_help .bubar.stopm "Stop mesh generation"

proc help_item { helptext } {p
    puts $helptext
}







proc show_help { } {
    
    set w .help
    
    if {[winfo exists .help] == 1} {
	wm withdraw $w
	wm deiconif $w
	focus $w 
    } {
	
	toplevel $w
	
	frame $w.buttons
	pack $w.buttons -side bottom -fill x -pady 2m
	button $w.buttons.done -text Done -command "destroy $w"
	pack $w.buttons.done  -side left -expand 1

	text $w.text -yscrollcommand "$w.scroll set" -setgrid true \
	    -width 60 -height 24 -wrap word
	scrollbar $w.scroll -command "$w.text yview"
	pack $w.scroll -side right -fill y
	pack $w.text -expand yes -fill both

    }
    $w.text configure -state normal
    $w.text delete 1.0 end
}






set bold "-background #43ce80 -relief raised -borderwidth 1"
set normal "-background {} -relief flat"


proc help_main { } {

    show_help;
    set w .help
    global bold
    global normal


    
    $w.text insert 0.0 \
	{NETGEN Help}
    $w.text insert end \n\n
    $w.text insert end \
	{1. General} d1
    $w.text insert end \n\n
    $w.text insert end \
	{2. Menu items } d2
    $w.text insert end \n\n

    foreach tag {d1 d2} {
	$w.text tag bind $tag <Any-Enter> "$w.text tag configure $tag $bold"
	$w.text tag bind $tag <Any-Leave> "$w.text tag configure $tag $normal"
    }
    
    $w.text tag bind d1 <1> { puts "general"; help_general }
    $w.text tag bind d2 <1> { help_menus }

    $w.text configure -state disabled
}






proc help_general { } {

    show_help;
    set w .help
    global bold
    global normal

    puts "general called"

    $w.text insert 0.0 \
	{NETGEN is an automatic three dimensional tetrahedral mesh generation system. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STEP or STL file format. NETGEN contains modules for mesh optimization and hierarchical mesh refinement.}

    $w.text configure -state disabled
}





proc help_menus { } {

    show_help;
    set w .help
    global bold
    global normal


    $w.text insert 0.0 \
	{The NETGEN Menu items are}
    $w.text insert end \n\n
    $w.text insert end \
	{1. File} d1
    $w.text insert end \n\n
    $w.text insert end \
	{2. Geometry } d2
    $w.text insert end \n\n
    $w.text insert end \
	{3. Mesh } d3
    $w.text insert end \n\n
    $w.text insert end \
	{4. View } d4
    $w.text insert end \n\n
    $w.text insert end \
	{5. Mesh-size } d5
    $w.text insert end \n\n
    $w.text insert end \
	{6. STL } d6

    foreach tag {d1 d2 d3 d4 d5 d6} {
	$w.text tag bind $tag <Any-Enter> "$w.text tag configure $tag $bold"
	$w.text tag bind $tag <Any-Leave> "$w.text tag configure $tag $normal"
    }
    
    $w.text tag bind d1 <1> {puts "File menu"}
    $w.text tag bind d2 <1> {puts "Geometry menu"}
    $w.text tag bind d3 <1> {puts "Mesh menu"}
    $w.text tag bind d4 <1> {puts "View menu"}
    $w.text tag bind d5 <1> {puts "Mesh-size menu"}
    $w.text tag bind d6 <1> {puts "STL menu"}

    $w.text configure -state disabled
}



