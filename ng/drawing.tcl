#
#  Creates a drawing frame, and binds mouse events
#
set oldmousex 0
set oldmousey 0
#


# if { 1 } {

# use this one for Togl 2.0
#if {[catch {togl .ndraw -width 400 -height 300  -rgba true -double true -depth true -privatecmap false -stereo false -indirect false }] } {

# if {[catch {togl .ndraw -width 400 -height 300  -rgba true -double true -depth true -privatecmap false -stereo false -indirect false -create init  -display draw -reshape reshape  }] } {    
#   puts "no OpenGL" 
# } {

set toglversion [Ng_GetToglVersion]
puts "togl-version : $toglversion"

set toglok 0
if { [Ng_GetToglVersion] == 2 } {
    # Togl 2.x
    if {[catch {togl .ndraw -width 400 -height 300  -rgba true -double true -depth true -privatecmap false -stereo false -indirect false  -create init  -display draw -reshape reshape }] } {    
        puts "no OpenGL" 
    } {
        # puts "have Togl 2.1"        
        set toglok 1
    }
} {
    # Togl 1.7
    if {[catch {togl .ndraw -width 400 -height 300  -rgba true -double true -depth true -privatecmap false -stereo false -indirect false  }] } {    
        puts "no OpenGL" 
    } {
        # puts "have Togl 1.7"
        set toglok 1
    }
}

if { $toglok == 1} {
    
    #
    pack .ndraw -expand true -fill both -padx 10 -pady 10
    catch { tkdnd::drop_target register .ndraw DND_Files }
    #
    bind .ndraw <Button-1> {
	set oldmousex %x; set oldmousey %y;
    }
    bind .ndraw <Button-2> {
	set oldmousex %x; set oldmousey %y;
    }
    bind .ndraw <Button-3> {
	set oldmousex %x; set oldmousey %y;
    }
    bind .ndraw <B1-Motion> {
	Ng_MouseMove $oldmousex $oldmousey %x %y $drawmode
	.ndraw render
	set oldmousex %x; set oldmousey %y;
    }

    bind .ndraw <Double-1> {
	Ng_MouseDblClick %x %y
	.ndraw render
	if { [winfo exists .bcprop_dlg] } { bcpropdialog }
	if { [winfo exists .surfacemeshsize_dlg] } { surfacemeshsizedialog }
	if { [winfo exists .fieldlines_dlg] } { fieldlinesdialog }
    }

    bind .ndraw <B2-Motion> {
	Ng_MouseMove $oldmousex $oldmousey %x %y move
	.ndraw render
	set oldmousex %x; set oldmousey %y;
    }

    bind .ndraw <B3-Motion> {
        if { $tcl_platform(os) == "Darwin" } {
            Ng_MouseMove $oldmousex $oldmousey %x %y move
        } { 
            Ng_MouseMove $oldmousex $oldmousey %x %y zoom
        }
	.ndraw render
	set oldmousex %x; set oldmousey %y;
    }
}


proc popupcheckredraw { vari { x 0 } } {
    upvar $vari varname
    if { $varname == 1 } {
	set varname 0
    } {
        #        puts "popup-redraw $vari"
	Ng_Vis_Set parameters
	redraw
    }
}
proc popupcheckredraw2 { vari boolvar { x 0 } } {
    upvar $vari varname
    if { $varname == 1 } {
	set varname 0
    } {
	Ng_SetVisParameters
	Ng_Vis_Set parameters
	if { $boolvar == 1 } { redraw }
	Ng_SetVisParameters
    }
}
proc popupcheckredraw3 { vari { x 0 } } {
    upvar $vari varname
    if { $varname == 1 } {
	set varname 0
    } {
	Ng_Vis_Set parameters
    }
}




proc redraw { {x 0} } {
    if {[winfo exists .ndraw]} { .ndraw render } 
}



bind . <Left> { Ng_MouseMove 0 0 -10 0 rotate; redraw }
bind . <Right> { Ng_MouseMove 0 0 10 0 rotate; redraw }
bind . <Up> { Ng_MouseMove 0 0 0 -10 rotate; redraw }
bind . <Down> { Ng_MouseMove 0 0 0 10 rotate; redraw }
bind . <Shift-Left> { Ng_MouseMove 0 0 -10 0 move; redraw }
bind . <Shift-Right> { Ng_MouseMove 0 0 10 0 move; redraw }
bind . <Shift-Up> { Ng_MouseMove 0 0 0 -10 move; redraw }
bind . <Shift-Down> { Ng_MouseMove 0 0 0 10 move; redraw }
bind . <Control-Up> { Ng_MouseMove 0 0 0 -10 zoom; redraw }
bind . <Control-Down> { Ng_MouseMove 0 0 0 10 zoom; redraw }

bind . <Button-4> \
   {event generate [focus -displayof %W] <MouseWheel> -delta  120}

bind . <Button-5> \
   {event generate [focus -displayof %W] <MouseWheel> -delta -120}

bind . <MouseWheel> { Ng_MouseMove 0 0 0 [expr {%D/-5}] zoom; redraw }



# Drop callbacks:
bind .ndraw <<Drop:DND_Files>> {
  #tk_messageBox -message "Dropped files: \"[join %D {, }]\""
  #%W state !active  
  set index [string first . %D end-6]
  #tk_messageBox -message $index
  set type [string range %D $index+1 end]
  # tk_messageBox -message $type
  set ispde [string match -nocase $type "pde"]
  set isgeo [expr max([string match -nocase $type "geo"],[string match -nocase $type "in2d"])]
  set ismesh [expr max([string match -nocase $type "vol"],[string match -nocase $type "vol.gz"])]
  set ispy [string match -nocase $type "py"]
  if {$ispde == 1} {
      AddRecentNGSFile %D;
      NGS_LoadPDE  %D;
      SetNumProcHelpMenu
      set selectvisual mesh;
      Ng_SetVisParameters
  }
  if {$ispy == 1} {
      AddRecentPYNGSFile %D;
      NGS_LoadPy  %D;
  }
  if {$isgeo == 1} {
	    AddRecentFile %D;
	    Ng_LoadGeometry %D;
	    Ng_ParseGeometry
	    set selectvisual geometry
	    Ng_SetVisParameters
	    redraw
	    wm title . [concat "$progname - " %D]
	    set dirname [file dirname %D]
	    set basefilename [file tail [file rootname %D]]
	}
	if {$ismesh == 1} {
	    AddRecentMeshFile %D;
	    Ng_LoadMesh %D; 
	    set selectvisual mesh
	    Ng_SetVisParameters
	    redraw
	    Ng_ReadStatus; 
	    wm title . [concat "$progname - " %D] 
	    set dirname [file dirname %D]
	    set basefilename [file tail [file rootname %D]]
	}
  return %A
}
