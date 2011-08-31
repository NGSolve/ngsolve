if { [catch { load libstlvis[info sharedlibextension] Ng_STL } result ] } {
#    puts "cannot load stl" 
#    puts "error: $result"
}


.ngmenu.geometry add separator

.ngmenu.geometry add command -label "STL Doctor..." \
    -command { stldoctordialog; }

.ngmenu.geometry add command -label "STL Info" \
    -command {
	set notriangles 0
	set minx 0
	set maxx 0
	set miny 0
	set maxy 0
	set minz 0
	set maxz 0
	set trigscons 0
	Ng_STLInfo notriangles minx maxx miny maxy minz maxz trigscons
	set msgtext "NO STL-Triangles : $notriangles\nGeometry:\nX = $minx - $maxx\nY = $miny - $maxy\nZ = $minz - $maxz\nConsistency Check = $trigscons\n"
	set msgtext "$msgtext Status: [Ng_STLInfo status]"
	tk_messageBox -title "STL Info" -message  $msgtext -type ok 
    }

