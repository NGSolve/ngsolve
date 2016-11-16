# netgen menus:

menu .ngmenu -tearoff 0  -relief raised -bd 2
. configure -menu .ngmenu

.ngmenu add cascade -label "File" -menu .ngmenu.file -underline 0
.ngmenu add cascade -label "Geometry" -menu .ngmenu.geometry -underline 0
.ngmenu add cascade -label "Mesh" -menu .ngmenu.mesh -underline 0
.ngmenu add cascade -label "View" -menu .ngmenu.view -underline 0
.ngmenu add cascade -label "Refinement" -menu .ngmenu.meshsize -underline 5

if { $userlevel == 3} {
    .ngmenu add cascade -label "Special" -menu .ngmenu.special -underline 3
}

.ngmenu add cascade -label "Help" -menu .ngmenu.help -underline 0


#####################################################
#                                                   #
#     Menu File                                     #
#                                                   #
#####################################################

menu .ngmenu.file

.ngmenu.file add command -label "Load Geometry..." -accelerator "<l><g>" \
    -command { 
	set types {
	    {"All Geometry types"   { .stl .stlb .step .stp .geo .in2d .igs .iges .brep .sat} }
	    {"IGES Geometry"	{.igs .iges} }
	    {"BREP OpenCascade Geometry"    {.brep} }
	    {"STL Geometry"        {.stl} }
	    {"Binary STL Geometry"    {.stlb} }
	    {"STEP Geometry"    {.step .stp} }
	    {"Geometry file"       {.geo} }
	    {"2D Geometry"   {.in2d } } 
	} 

	set ACISavailable [Ng_ACISCommand isACISavailable]
	if {$ACISavailable == "yes" } {
	    lappend types {"ACIS Geometry" {.sat} }
	}

	if {[catch {
	    set file [tk_getOpenFile -filetypes $types -initialdir $dirname -typevariable loadgeomtypevar]
	}]} {
	    set file [tk_getOpenFile -filetypes $types -initialdir $dirname]
	}

	if {$file != ""} {
	    AddRecentFile $file
	    Ng_LoadGeometry $file 
	    Ng_ParseGeometry
#	    if { [Ng_STLInfo status]=="ERROR" } {
#		tk_messageBox -message  "STL ERROR: \n [Ng_STLInfo statustext]" -type ok
#	    }
	    set selectvisual geometry
	    Ng_SetVisParameters
	    redraw
	    wm title . [concat "$progname - " $file]
	    set dirname [file dirname $file]
	    set basefilename [file tail [file rootname $file]]

	    if { $hasocc == "yes" } {
		rebuildoccdialog
	    }
	}
    }



.ngmenu.file add command -label "Save Geometry..." \
    -command { 
	set occgeometryloaded [Ng_OCCCommand isoccgeometryloaded]
	puts $occgeometryloaded
	if {$occgeometryloaded == 1 } {
	    set types {
		{"IGES Geometry file"   {.igs} } 
		{"STEP Geometry file"   {.stp} } 
		{"STL Geometry file"   {.stl} } 
		{"STL BIN Geometry file"   {.stlb} } 
	    }
	} {
	    set types {
		{"STL Geometry file"   {.stl} } 
		{"STL BIN Geometry file"   {.stlb} } 
	    }
	}

	set ACISavailable [Ng_ACISCommand isACISavailable]
	puts $ACISavailable
	if {$ACISavailable == "yes" } {
	    lappend types {"ACIS Geometry" {.sat} }
	}

	set file [tk_getSaveFile -filetypes $types -initialdir $dirname -initialfile $basefilename ]
	if {$file != ""} {
	    Ng_SaveGeometry $file 
	}
    }
 


.ngmenu.file add cascade -label "Recent Files" -menu .ngmenu.file.recent 
menu .ngmenu.file.recent -tearoff 0


proc AddRecentFile { filename } {
    global progname
    global dirname
    catch { [.ngmenu.file.recent delete $filename] }
    .ngmenu.file.recent insert 0 command -label $filename \
	-command "AddRecentFile {$filename}; 
                  Ng_LoadGeometry {$filename}; 
		  Ng_ParseGeometry;
		  set selectvisual geometry;
		  Ng_SetVisParameters;
	          redraw;
		  wm title . [concat \" $progname - $filename \"];
                  set dirname {[file dirname $filename]};
                  set basefilename {[file tail [file rootname $filename]]};
        	  rebuildoccdialog;"

    
    if { [.ngmenu.file.recent index last] >= 6 } {
	.ngmenu.file.recent delete last }
    
    saveinifile;
    }
loadinifile;

.ngmenu.file add separator


.ngmenu.file add command -label "Load Mesh..." -accelerator "<l><m>" \
    -command {
	set types {
	    {"Mesh file"   {.vol .vol.gz}	} }
	set file [tk_getOpenFile -filetypes $types -defaultextension ".vol"]
	if {$file != ""} {
	    AddRecentMeshFile $file;
	    Ng_LoadMesh $file; 
	    set selectvisual mesh
	    Ng_SetVisParameters
	    redraw
	    Ng_ReadStatus; 
#	    Ng_MeshSizeFromSurfaceMesh
	    wm title . [concat "$progname - " $file] 
	    set dirname [file dirname $file]
	    set basefilename [file tail [file rootname $file]]
	}
    }



# astrid
.ngmenu.file add cascade -label "Recent Meshes" -menu .ngmenu.file.recentmesh 
menu .ngmenu.file.recentmesh


proc AddRecentMeshFile { filename } {
    global progname
    global dirname
    catch { [.ngmenu.file.recentmesh delete $filename] }
    .ngmenu.file.recentmesh insert 0 command -label $filename \
	-command "AddRecentMeshFile {$filename}; 
                  Ng_LoadMesh {$filename};
		  set selectvisual mesh;
		  Ng_SetVisParameters;
	          redraw;
		  wm title . [concat \" $progname - $filename \"];
                  set dirname {[file dirname $filename]};
                  set basefilename {[file tail [file rootname $filename]]};
                  rebuildoccdialog;"
    
    if { [.ngmenu.file.recentmesh index last] >= 6 } {
	.ngmenu.file.recentmesh delete last }
   
    savemeshinifile;
    }
loadmeshinifile;

# astrid ende


.ngmenu.file add command -label "Save Mesh..." -accelerator "<s><m>" \
    -command {
	set types {
	    {"Mesh file"   {.vol .vol.gz}	} }

	set file [tk_getSaveFile -filetypes $types -defaultextension ".vol.gz" -initialfile $basefilename -initialdir $dirname ]
	if {$file != ""} {
	    Ng_SaveMesh $file }
	AddRecentMeshFile $file;

    }

.ngmenu.file add command -label "Merge Mesh..." \
    -command {
	set types {
	    {"Mesh file"   {.vol}	} }
	set file [tk_getOpenFile -filetypes $types -defaultextension ".vol"]
	if {$file != ""} {
	    Ng_MergeMesh $file; 
	    set selectvisual mesh
	    Ng_SetVisParameters
	    redraw
	    Ng_ReadStatus; 
	}
    }





.ngmenu.file add command -label "Import Mesh..." \
    -command { 
	set types {
	    {"Neutral format"  {.mesh .emt} }
	    {"Surface mesh format"  {.surf} }
	    {"Universal format"  {.unv} }
	    {"Olaf format"  {.emt} }
	    {"TET format" {.tet} }
	    {"Pro/ENGINEER neutral format" {.fnf} }
	          }
	set file [tk_getOpenFile -filetypes $types ]
	if {$file != ""} {
	    Ng_ImportMesh $file 
	    set selectvisual mesh
	    Ng_SetVisParameters
	    redraw
	    Ng_ReadStatus; 
	}
    }


.ngmenu.file add command -label "Export Mesh..." \
    -command {

# 	global meshexportformats
	foreach exportformat $meshexportformats {
	    if { [lindex $exportformat 0] == $exportfiletype } {
		set extension [lindex $exportformat 1]
	    }
	}

	if { $exportfiletype == "Elmer Format"} {
	    set file [file nativename [tk_chooseDirectory -title "Elmer Mesh Export - Select Directory"]]
        } elseif { $exportfiletype == "OpenFOAM 1.5+ Format"} {
	    set file [file nativename [tk_chooseDirectory -title "OpenFOAM 1.5+ Mesh Export - Select Case Directory"]]
        } elseif { $exportfiletype == "OpenFOAM 1.5+ Compressed"} {
	    set file [file nativename [tk_chooseDirectory -title "OpenFOAM 1.5+ Mesh Export - Select Case Directory"]]
        } else {
#	    set file [tk_getSaveFile  -filetypes "{ \"$exportfiletype\" {$extension} }" ]
	    set file [tk_getSaveFile  -filetypes "{ \"$exportfiletype\" {*}}" ]
	}

	if {$file != ""} {
	    Ng_ExportMesh $file $exportfiletype 
	}
    }

.ngmenu.file add cascade -label "Export Filetype" -menu .ngmenu.file.filetype 

menu .ngmenu.file.filetype 


.ngmenu.file add separator


.ngmenu.file add command -label "Save Solution..." \
    -command { 
	set types { 
            {"Solution File"  {.sol} } 
            {"VTK File"  {.vtk} } 
        }
	set file [tk_getSaveFile -filetypes $types ]
	if {$file != ""} {
	    Ng_SaveSolution $file 
	}
    }
#-defaultextension ".sol"  ]

.ngmenu.file add command -label "Import Solution..." \
    -command { 
	set types { {"Solution File"  {.sol} } }
	set file [tk_getOpenFile -filetypes $types -defaultextension ".sol"  ]
	if {$file != ""} {
	    Ng_ImportSolution $file 
	    set selectvisual solution
	    Ng_SetVisParameters
	    redraw
	}
    }






set demostarttime [clock clicks -millisecond]
set stopdemo 0
proc demoredraw { } {
    global demostarttime
    global stopdemo
    set curtime [clock clicks -millisecond]
    set result [ Ng_DemoSetTime [expr $curtime - $demostarttime] ]
    redraw
    global videoactive
    if { $videoactive == 1 } {
        puts "addframe"
        .ndraw Ng_VideoClip addframe
    }
    if { $result == 0 && $stopdemo == 0 } {
	after 1 { demoredraw }
    }
}
.ngmenu.file add command -label "Show Demo..." \
    -command {
	set types { {"Demo File"  {.dem} } }
	set file [tk_getOpenFile -filetypes $types -defaultextension ".dem"  ]
	if {$file != ""} {
	    Ng_ShowDemo $file 
	    set demostarttime [clock clicks -millisecond]
	    set stopdemo 0
	    demoredraw
 	}
     }




.ngmenu.file add separator

.ngmenu.file add command -label "Snapshot..." \
    -command { 
	set types { 
	    {"JPG file" {.jpg} } 
	    {"GIF file" {.gif} } 
	    {"PPM file" {.ppm} } 
	}
	set file [tk_getSaveFile -filetypes $types]
#  -defaultextension ".ppm"]
	if {$file != ""} {
	    .ndraw Ng_SnapShot $file }
    }


.ngmenu.file add cascade -label "Video clip" -menu .ngmenu.file.video
menu .ngmenu.file.video

set videoactive 0
.ngmenu.file.video add command -label "start..." \
    -command { 
 	set types { 
 	    {"MPG file" {.mpg} } 
 	}
 	set file [tk_getSaveFile -filetypes $types]
 	if {$file != ""} {
 	    .ndraw Ng_VideoClip init $file 
            global videoactive
            set videoactive 1
        }
     }

.ngmenu.file.video add command -label "add frame..." \
    -command {.ndraw Ng_VideoClip addframe }

.ngmenu.file.video add command -label "one cycle" \
    -command {
	set visoptions.redrawperiodic 1
	for { set j 0 } { $j < 100 } { incr j } {
	    puts "j =  $j"
	    Ng_Vis_Set time [expr (1000 * $j / 100)]
	    redraw
	    .ndraw Ng_VideoClip addframe 
	    after 200
	}
    }

.ngmenu.file.video add command -label "finalize..." \
    -command {
        .ndraw Ng_VideoClip finalize 
        global videoactive
        set videoactive 0
    }



.ngmenu.file add command -label "Save Options" \
    -command { saveoptions }


    

.ngmenu.file add separator


## herbert tcl load menue
# .ngmenu.file add command -label "Run tests ..." \
\#    -command { runtestdialog }
##
# .ngmenu.file add separator

.ngmenu.file add command -label "Quit" -accelerator "<q>" \
    -command { 
        puts "Thank you for using $progname"; 

        if { [catch { unload libngsolve[info sharedlibextension] ngsolve } result ] } {
            # puts "cannot unload ngsolve" 
            # puts "error: $result"
        } 

        Ng_Exit; 
        destroy . 
    }
# exit


#####################################################
#                                                   #
#     Menu Mesh                                     #
#                                                   #
#####################################################

menu .ngmenu.mesh
.ngmenu.mesh add command -label "Generate Mesh" -accelerator "<g><m>" \
    -command { 
	set selectvisual mesh
	Ng_SetVisParameters
	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}
	Ng_ReadStatus
	redraw
    }

.ngmenu.mesh add command -label "Stop Meshing" \
    -command { Ng_StopMeshing }

.ngmenu.mesh add command -label "Meshing Options..." \
    -command meshingoptionsdialog


.ngmenu.mesh add separator

.ngmenu.mesh add command -label "Delete Mesh" \
    -command { Ng_New mesh; Ng_ReadStatus; redraw }

.ngmenu.mesh add command -label "Delete Vol Mesh" \
    -command { Ng_DeleteVolMesh; Ng_ReadStatus; redraw }


.ngmenu.mesh add command -label "Mesh Info" \
    -command {
	set dim [Ng_MeshInfo dim]
	set np [Ng_MeshInfo np]
	set ne [Ng_MeshInfo ne]
	set nse [Ng_MeshInfo nse]
	set nseg [Ng_MeshInfo nseg]
	set bbox [Ng_MeshInfo bbox]
	tk_messageBox -message  "Dimension: $dim\nPoints: $np\nElements: $ne\nSurface Els: $nse\nSegments: $nseg\nxmin [lindex $bbox 0] xmax [lindex $bbox 1]\nymin [lindex $bbox 2] ymax [lindex $bbox 3]\nzmin [lindex $bbox 4] zmax [lindex $bbox 5]"
    }


.ngmenu.mesh add command -label "Mesh Quality" \
    -command {
	set inplanemin 0
	set inplanemax 0
	set betplanemin 0
	set betplanemax 0
	Ng_MeshQuality inplanemin inplanemax betplanemin betplanemax
	puts "Triangle angles : $inplanemin - $inplanemax"
	puts "Tet angles      : $betplanemin - $betplanemax"
	tk_messageBox -message  "Triangle angles : $inplanemin - $inplanemax \n Tet angles      : $betplanemin - $betplanemax"
    }

# .ngmenu.mesh add command -label "Quality Plot" \
#    -command { qualityviewdialog 1 }




.ngmenu.mesh add command -label "Check Surface Mesh" \
    -command { Ng_CheckSurfaceMesh }
.ngmenu.mesh add command -label "Check Volume Mesh" \
    -command { Ng_CheckVolumeMesh }

.ngmenu.mesh add command -label "Edit Boundary Conditions..." \
    -command { bcpropdialog }

if { $userlevel == 3 } {
    .ngmenu.mesh add command -label "Mesh Doctor..." \
	-command { meshdoctordialog }
}

.ngmenu.mesh add command -label "METIS Mesh Partitioning..." \
	-command { METISdialog }

.ngmenu.mesh add separator

.ngmenu.mesh add command -label "Analyze Geometry" \
    -command { Ng_GenerateMesh ag ag; Ng_ReadStatus; redraw }
.ngmenu.mesh add command -label "Mesh Edges" \
    -command { Ng_GenerateMesh me me; Ng_ReadStatus; redraw }
.ngmenu.mesh add command -label "Mesh Surface" \
    -command { set selectvisual mesh; Ng_SetVisParameters; \
		   Ng_GenerateMesh ms ms; Ng_ReadStatus; redraw }

.ngmenu.mesh add command -label "Optimize Surface" \
    -command { Ng_GenerateMesh os os cmsmSm; redraw }

.ngmenu.mesh add cascade -label "Surface Optim. Step" -menu .ngmenu.mesh.surfoptstep 

menu .ngmenu.mesh.surfoptstep 
.ngmenu.mesh.surfoptstep add command -label "Mesh Smoothing" \
    -command { Ng_GenerateMesh os os m; redraw}
.ngmenu.mesh.surfoptstep add command -label "Edge swapping (topologic)" \
    -command { Ng_GenerateMesh os os s; redraw}
.ngmenu.mesh.surfoptstep add command -label "Edge swapping (metric)" \
    -command { Ng_GenerateMesh os os S; redraw}
.ngmenu.mesh.surfoptstep add command -label "Combine points" \
    -command { Ng_GenerateMesh os os c; redraw}


.ngmenu.mesh add separator
.ngmenu.mesh add command -label "Mesh Volume" \
    -command { Ng_GenerateMesh mv mv; Ng_ReadStatus }
.ngmenu.mesh add command -label "Optimize Volume" \
    -command { Ng_GenerateMesh ov ov; Ng_ReadStatus }
.ngmenu.mesh add command -label "Smooth Opt Volume" \
    -command { Ng_GenerateMesh ov ov m; Ng_ReadStatus }
.ngmenu.mesh add command -label "Smooth Opt Volume Jacobian" \
    -command { Ng_GenerateMesh ov ov j; Ng_ReadStatus }



#####################################################
#                                                   #
#     Menu Geometry                                 #
#                                                   #
#####################################################

menu .ngmenu.geometry







#####################################################
#                                                   #
#     Menu View                                     #
#                                                   #
#####################################################

menu .ngmenu.view
.ngmenu.view add command -label "Zoom all" \
    -command { Ng_ZoomAll; redraw }
.ngmenu.view add command -label "Center" \
    -command { Ng_Center; redraw }

.ngmenu.view add command -label "x-y plane" \
    -command { Ng_StandardRotation xy; redraw }
.ngmenu.view add command -label "y-x plane" \
    -command { Ng_StandardRotation yx; redraw }
.ngmenu.view add command -label "x-z plane" \
    -command { Ng_StandardRotation xz; redraw }
.ngmenu.view add command -label "z-x plane" \
    -command { Ng_StandardRotation zx; redraw }
.ngmenu.view add command -label "y-z plane" \
    -command { Ng_StandardRotation yz; redraw }
.ngmenu.view add command -label "z-y plane" \
    -command { Ng_StandardRotation zy; redraw }

.ngmenu.view add command -label "Viewing Options..." \
    -command { viewingoptionsdialog; redraw }
.ngmenu.view add command -label "Clipping Plane..." \
    -command { clippingdialog; redraw }
.ngmenu.view add command -label "Solution Data..." \
    -command { visual_dialog; redraw }
.ngmenu.view add checkbutton -variable viewqualityplot \
    -label "Quality Plot" \
    -command { qualityviewdialog $viewqualityplot }
.ngmenu.view add checkbutton -variable memuseplot \
    -label "Memory Usage" \
    -command { memusedialog $memuseplot }




#####################################################
#                                                   #
#     Menu Refinement                               #
#                                                   #
#####################################################
#
# Mesh size menu
#
menu .ngmenu.meshsize
.ngmenu.meshsize add command -label "Refine uniform" \
    -command { Ng_Refine; Ng_HighOrder ${options.elementorder}; Ng_ReadStatus; redraw }

.ngmenu.meshsize add command -label "Second Order" \
    -command { Ng_SecondOrder; Ng_ReadStatus; redraw }

.ngmenu.meshsize add command -label "Validate Second Order" \
    -command { Ng_ValidateSecondOrder; Ng_ReadStatus; redraw }

.ngmenu.meshsize add command -label "High Order" \
    -command { Ng_HighOrder ${options.elementorder}; Ng_ReadStatus; redraw }

.ngmenu.meshsize add separator

.ngmenu.meshsize add command -label "Refinement Dialog..." \
    -command { refinementdialog }
.ngmenu.meshsize add command -label "Load Meshsize..." \
    -command {
	set types {
	    {"Meshsize file"   {.msz}	} }
	set file [tk_getOpenFile -filetypes $types]
	if {$file != ""} {
	    Ng_LoadMeshSize $file; 
	}
    }
.ngmenu.meshsize add command -label "MS from Surf Mesh" \
    -command { Ng_MeshSizeFromSurfaceMesh }


if { $userlevel == 3 } {
.ngmenu.meshsize add command -label "Singular point ms" \
    -command { Ng_SingularPointMS; }

.ngmenu.meshsize add command -label "Singular edge ms" \
    -command { Ng_SingularEdgeMS; }

.ngmenu.meshsize add separator

set bisectfilename "";

.ngmenu.meshsize add command -label "Bisection" \
    -command { Ng_ReadStatus; set oldnp 0; set newnp $status_np; 
#	Ng_BisectCopyMesh; 
#	Ng_Split2Tets;
	Ng_ReadStatus;
	
	while { $oldnp < $newnp } {
#	    if { $level == 0 } {
#		Ng_ExportMesh feppmesh.vol fepp;
#	    } {
#		Ng_ExportMesh feppmesh$level feppml 
#	    }
	    set level [expr $level+1]
	    if { $bisectfilename == ""} {
		Ng_Bisect;
	    } else {
		Ng_Bisect $bisectfilename;
	    }
#	    Ng_HighOrder ${options.elementorder} "noparallel"
#	    Ng_Split2Tets;
	    Ng_ReadStatus;
	    redraw; 
	    
	    if { $bisectfilename == ""} {
		set oldnp $newnp;
		set newnp $status_np;
		puts "oldnp $oldnp newnp $newnp";
	    } else {
		set oldnp $newnp;
	    }
	}
    }
#    -command { Ng_Bisect; Ng_ReadStatus; redraw }
#    -command { exec netgen abc >outfile 2>errfile; Ng_ReadStatus; redraw }

}

.ngmenu.meshsize add command -label "Load Refinement Info..." \
    -command {
	set types {
	    {"Refinement info" {.refine} }}
	set bisectfilename [tk_getOpenFile -filetypes $types]
    }

.ngmenu.meshsize add command -label "Z-Refinement" \
    -command { Ng_ZRefinement 2; Ng_ReadStatus; redraw }


# .ngmenu.meshsize add command -label "hp-Refinement" \
\#    -command { Ng_HPRefinement 4; Ng_ReadStatus; redraw }

.ngmenu.meshsize add cascade -label "hp-Refinement" -menu .ngmenu.meshsize.hpref
menu .ngmenu.meshsize.hpref
.ngmenu.meshsize.hpref add command -label "1 Level" \
    -command { Ng_HPRefinement 1; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "2 Levels" \
    -command { Ng_HPRefinement 2; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "3 Levels" \
    -command { Ng_HPRefinement 3; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "4 Levels" \
    -command { Ng_HPRefinement 4; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "5 Levels" \
    -command { Ng_HPRefinement 5; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "6 Levels" \
    -command { Ng_HPRefinement 6; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "7 Levels" \
    -command { Ng_HPRefinement 7; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "8 Levels" \
    -command { Ng_HPRefinement 8; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "9 Levels" \
    -command { Ng_HPRefinement 9; Ng_ReadStatus; redraw }
.ngmenu.meshsize.hpref add command -label "10 Levels" \
    -command { Ng_HPRefinement 10; Ng_ReadStatus; redraw }


.ngmenu.meshsize add command -label "Split to Tets" \
    -command { Ng_Split2Tets; Ng_ReadStatus; redraw }





#####################################################
#                                                   #
#     Menu Special                                  #
#                                                   #
#####################################################

menu .ngmenu.special
.ngmenu.special add command -label "Prismatic Boundary Layer" \
    -command { Ng_GenerateBoundaryLayer; redraw }
.ngmenu.special add command -label "Insert virtual boundary layer" \
    -command { Ng_InsertVirtualBL; redraw }
.ngmenu.special add command -label "Cut off and combine with other" \
    -command { 
	set types { {"Mesh file"   {.vol}	} }
	set file [tk_getOpenFile -filetypes $types]
	if {$file != ""} {
	    Ng_CutOffAndCombine $file;  }
	redraw 
    }
.ngmenu.special add command -label "Helmholtz Mesh grading" \
    -command { Ng_HelmholtzMesh; }
.ngmenu.special add cascade -label "Colour-based boundary conditions" -menu .ngmenu.special.colbndcond

menu .ngmenu.special.colbndcond 
 .ngmenu.special.colbndcond add command -label "Inspect Colours in mesh" \
    -command { currmeshcoloursdialog }
    
 .ngmenu.special.colbndcond add separator	
    
 .ngmenu.special.colbndcond add command -label "Automatic Assignment" \
    -command { Ng_AutoColourBcProps auto; redraw }
	
 .ngmenu.special.colbndcond add separator	

 set ocffile [file join ${ngdir} netgen.ocf];
 
 .ngmenu.special.colbndcond add command -label "Select Colour Profile file" \
    -command {
	set types { {"Colour Profile file"   {.ocf}   } }
	set ocffile [tk_getOpenFile -filetypes $types]
	if {$ocffile == ""} {
	    set ocffile [file join ${ngdir} netgen.ocf]; }
	} 
 .ngmenu.special.colbndcond add command -label "Profile based Assignment" \
	-command { Ng_AutoColourBcProps profile ${ocffile}; redraw }


# menu .mbar.stl.menu
# .mbar.stl.menu add command -label "STL options" \
#     -command { stloptionsdialog; }
#.mbar.stl.menu add command -label "STL Doctor" \
#    -command { stldoctordialog; }


#####################################################
#                                                   #
#     Menu Help                                     #
#                                                   #
#####################################################




menu .ngmenu.help
# .ngmenu.help add command -label "Ng Help..." \
\#	-command { help_main }
# .ngmenu.view add checkbutton -variable showsensitivehelp \
#	-label "Sensitve Help" \
#	-command { sensitivehelpdialog $showsensitivehelp }
.ngmenu.view add checkbutton -label "Help Line" -variable showhelpline \
	-command {
    if { $showhelpline == 1} {
	pack .helpline -before .statbar -side bottom -fill x -padx 3p
    } {
	pack forget .helpline 
    }
} 

.ngmenu.help add command -label "About..." \
    -command {
tk_messageBox -message "This is NETGEN \nmainly written by \nJoachim Schoeberl \nthanks to \nRobert Gaisbauer, Johannes Gerstmayr, Philippose Rajan"
}

# tk_menuBar .mbar .mbar.file .mbar.mesh .mbar.test .mbar.help
# focus .mbar




#####################################################
#                                                   #
#     Button bar                                    #
#                                                   #
#####################################################

ttk::frame .bubar
# -relief raised
#-relief raised -bd 2
pack .bubar -side top -fill x

button .bubar.testb -text "Test" -command { Ng_SaveGeometry }
ttk::button .bubar.surfm -text "Generate Mesh" -command \
    { 
	.ngmenu.mesh invoke "Generate Mesh"; 
#	set selectvisual mesh; 
#	Ng_SetVisParameters;
#	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}
#	redraw 
    }
ttk::button .bubar.stopm -text "Stop" -command \
    { 
	# Ng_StopMeshing;  
	set multithread_terminate 1;
	set stopdemo 1;
    }
ttk::button .bubar.exitb -text "Quit" \
    -command { 
	   set ans [tk_messageBox -title "Quit Netgen?" -message "Do you really want to quit Netgen?" -type yesno -default "no" -icon question]
	   if { $ans == "yes" } {
	     .ngmenu.file invoke "Quit"; 
	   }	 
	}
pack  .bubar.exitb .bubar.surfm .bubar.stopm -side left

#button .bubar.scan -text "Scan" \
#    -command { Ng_ParseGeometry; set selectvisual geometry; Ng_SetVisParameters; redraw }

ttk::button .bubar.zoomall -text "Zoom All" \
    -command { Ng_ZoomAll; redraw }

ttk::button .bubar.center -text "Center" \
    -command { Ng_Center; redraw }

# tk_optionMenu .bubar.modesel drawmode "rotate" "move  " "zoom  "
# tixOptionMenu .bubar.modesel \
     # -options {
 	# label.width  0
 	# label.anchor e
 	# menubutton.width 6
     # } \
     # -variable drawmode

# .bubar.modesel add command rotate -label Rotate
# .bubar.modesel add command move -label Move
# .bubar.modesel add command zoom -label Zoom

ttk::menubutton .bubar.modesel -menu .bubar.modesel.menu -text "" -width 6

menu .bubar.modesel.menu  -tearoff 0

.bubar.modesel.menu add command -label "Rotate" -command "set drawmode \"rotate\" ;.bubar.modesel configure -text \"Rotate\""
.bubar.modesel.menu add command -label "Move" -command "set drawmode \"move\" ;.bubar.modesel configure -text \"Move\""
.bubar.modesel.menu add command -label "Zoom" -command "set drawmode \"zoom\" ;.bubar.modesel configure -text \"Zoom\""

.bubar.modesel.menu invoke "Rotate"



      
set viewvals { geometry specpoints mesh solution}
if { $userlevel == 3} {
    set viewvals { geometry mesh specpoints surfmeshing modelview solution}
}

set viewvallabs(cross)     "Cross" 
set viewvallabs(geometry)  "Geometry" 
set viewvallabs(mesh)      "Mesh" 
set viewvallabs(specpoints) "Edges" 
set viewvallabs(surfmeshing) "Mesh Gen" 
set viewvallabs(modelview)     "Modeller" 
set viewvallabs(solution)     "Solution" 

# tixOptionMenu .bubar.selview \
    # -options {
	# label.width  0
	# label.anchor e
	# menubutton.width 10
    # } \

# foreach viewv $viewvals {
    # .bubar.selview add command $viewv -label $viewvallabs($viewv)
# }

# .bubar.selview config -variable selectvisual
# .bubar.selview config -command { Ng_SetVisParameters; redraw }



# pack .bubar.modesel -side right
# pack forget .bubar.modesel

pack .bubar.center .bubar.zoomall -side right
# pack .bubar.selview -side right

.ngmenu.view add checkbutton -variable viewrotatebutton \
    -label "Enable LeftButton Selection" \
    -command { 
	if { $viewrotatebutton } {
	    pack .bubar.modesel -side right
	} {
	    pack forget .bubar.modesel
	}
    }

menu .bubar.selviewmenu
ttk::menubutton .bubar.selview1 -menu .bubar.selviewmenu -text "Geometry"
foreach viewv $viewvals {
    .bubar.selviewmenu add command -label $viewvallabs($viewv) -command \
        ".bubar.selview1 configure -text \"$viewvallabs($viewv)\" ; set selectvisual $viewv ; Ng_SetVisParameters; redraw"
}
pack .bubar.selview1 -side right
# .bubar.selviewmenu invoke $viewvallabs($selectvisual)

trace add variable selectvisual write selvis_monitor
proc selvis_monitor { name args } {
    global selectvisual viewvallabs    
    .bubar.selviewmenu invoke $viewvallabs($selectvisual)
} 
# set selectvisual solution


#####################################################
#                                                   #
#     Status bar                                    #
#                                                   #
#####################################################

label .helpline -text "None"
pack forget .helpline -side bottom -fill x

ttk::frame .statbar -relief flat
# -bd 2
pack .statbar -side bottom -fill x

ttk::label .statbar.ptslabel -text "   Points: "
ttk::label .statbar.ptsval -textvariable status_np
ttk::label .statbar.elslabel -text "   Elements: "
ttk::label .statbar.elsval -textvariable status_ne
ttk::label .statbar.selslabel -text "   Surf Elements: "
ttk::label .statbar.selsval -textvariable status_nse
# label .statbar.memlabel -text "   Mem: "
# label .statbar.memval -textvariable mem_moveable
ttk::label .statbar.task -textvariable status_task

pack .statbar.ptslabel .statbar.ptsval -side left -ipady 3p 
pack .statbar.elslabel .statbar.elsval -side left -ipady 3p 
pack .statbar.selslabel .statbar.selsval -side left -ipady 3p

# if { $userlevel == 3 } {
#    pack .statbar.memlabel .statbar.memval -side left -ipady 3p
# }


#tixMeter .statbar.per -value 0 -text 0%
ttk::progressbar .statbar.per -value 0 -maximum 1
#.statbar.per configure -fillcolor blue

pack .statbar.per -side right
pack .statbar.task -side right -ipady 4

set qualbaraxis(0) 0
set qualbar(0) 0
set qualbarnull(0) 0



proc timer2 { } {
    global status_np
    global status_ne
    global status_nse
    global multithread_running
    global multithread_redraw
    global status_working
    global status_task
    global status_percent
    global status_tetqualclasses
    

    Ng_ReadStatus 

    if { $multithread_redraw == 1 } {
        # non-blocking redraw
	set multithread_redraw 0;
	redraw;
        
        global videoactive
        if { $videoactive == 1 } {
            puts "addframe"
            .ndraw Ng_VideoClip addframe
        }
    }
    if { $multithread_redraw == 2 } {
        # blocking redraw
	redraw;
	set multithread_redraw 0;
        
        global videoactive
        if { $videoactive == 1 } {
            puts "addframe"
            .ndraw Ng_VideoClip addframe
        }
        after 1 { timer2 }
        return
    }


    # global mem_moveable
    # set mem_moveable [Ng_MemInfo moveable]


    .statbar.per configure -value [expr $status_percent/100]
#    -text [format %2.1f [expr 0.1*int(10*$status_percent)]]%


    if { $multithread_running } {
	pack .statbar.per -side right -before .statbar.task -padx 6
    } { 
	pack forget .statbar.per
    }
	


    # tet quality
    if {[winfo exists .qualityview_dlg] == 1} {
	
	global qualbar
	global qualbarnull
	global qualbaraxis

	set maxval 0
	for {set i 0} {$i < 20} {incr i} {
	    if {[lindex $status_tetqualclasses $i] > $maxval} {
		set maxval [lindex $status_tetqualclasses $i]
	    }
	} 

	set ubound 1
	while { $ubound < $maxval } {
	    set ubound [expr {10 * $ubound}]
	}
	if { $ubound/5 > $maxval } {
	    set ubound [expr $ubound/5]
	}
	if { $ubound/2 > $maxval } {
	    set ubound [expr $ubound/2]
	}


	
	for {set i 1} {$i <= 5} {incr i} {
	    # global qualbaraxis($i)

	    set value [expr { $i * $ubound / 5 }]
	    .qualityview_dlg.c dchars $qualbaraxis($i) 0 end
	    .qualityview_dlg.c insert $qualbaraxis($i) end $value  
	}

	
	for {set i 0} {$i < 20} {incr i} {
	    set x1 [expr {100 + ($i*15) + 2}]
	    set x2 [expr {$x1+10}]
	    
	    set nbrs [lindex $status_tetqualclasses $i]
	    set y [expr (249 - (200 * $nbrs / $ubound ) )]
	    
	    # global qualbar($i)
	    .qualityview_dlg.c coords $qualbar($i) $x1 250 $x2 $y

#	    global qualbarnull($i)
	    if { $nbrs == 0 } {
		.qualityview_dlg.c itemconfigure $qualbarnull($i) -text 0
	    } {
		.qualityview_dlg.c itemconfigure $qualbarnull($i) -text "" 
	    }		
	}
	
    }


    if {[winfo exists .memuse_dlg] == 1} {    
	
	global memmark
	set usemb [Ng_MemInfo usedmb]
	for {set i 0} {$i < [string length $usemb] } {incr i} {
	    if { [string index $usemb $i] == 0 } {
		.memuse_dlg.c coords $memmark($i)  [expr 50+$i] 68 [expr 50+$i] 70
	    } {
		.memuse_dlg.c coords $memmark($i)  [expr 50+$i] 50 [expr 50+$i] 70
	    }
	}

    }
    after 30 { timer2 }
}
# after 1000 { timer2 }
timer2




proc bgerror { error } {
    global errorInfo userlevel
    if { $userlevel == 3} {
	puts "ERROR: $error" 
	puts "errinfo: $errorInfo"
    }
    tk_messageBox -title "Error Message" -message $error -type ok 
}






proc smh2 { menuitem } {
    if {[catch {$menuitem entrycget active -label} name]} {
	set name "    "
    } 
    show_menu_help $name 
    update idletasks
}

bind .ngmenu <<MenuSelect>> { smh2 %W }
bind .ngmenu.file <<MenuSelect>> { smh2 %W }
bind .ngmenu.geometry <<MenuSelect>> { smh2 %W }
bind .ngmenu.mesh <<MenuSelect>> { smh2 %W }
bind .ngmenu.view <<MenuSelect>> { smh2 %W }
bind .ngmenu.meshsize <<MenuSelect>> { smh2 %W }
bind .ngmenu.special <<MenuSelect>> { smh2 %W }
bind .ngmenu.help <<MenuSelect>> { smh2 %W }


# command bindings  
bind . <q> { .ngmenu.file invoke "Quit" }
bind . <l><g> { .ngmenu.file invoke "Load Geometry..." }  ; 
bind . <l><m> { .ngmenu.file invoke "Load Mesh..." }  ;
bind . <s><m> { .ngmenu.file invoke "Save Mesh..." }  ;
bind . <r><f> { .ngmenu.file activate "Recent Files" }  ;
bind . <n><p> { newprimitivedialog }      ; # 
bind . <e><p> { editprimitivedialog }
bind . <e><s> { newsoliddialog }
bind . <g><m> { .ngmenu.mesh invoke "Generate Mesh" }  ;





