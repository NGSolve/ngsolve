puts "Hello from Tcl"

set w .tcldemo;
toplevel $w
wm withdraw $w
wm geom $w +100+100
wm deiconify $w;
wm title $w "Tcl-Demo";

set evalfunction "sin(3*pi*x)"
tixLabelEntry $w.function -label "Enter function: "  \
    -labelside top \
    -options {  
        entry.textVariable evalfunction
        entry.width 35
        label.width 25
        label.anchor w
    } 

button $w.tcldemobutton -text "Draw" -command {
    puts "Tcl: function is $evalfunction";
    NGS_TclDemo $evalfunction
    Ng_Vis_Set parameters; redraw;
}


pack $w.function $w.tcldemobutton

focus .tcldemo;
