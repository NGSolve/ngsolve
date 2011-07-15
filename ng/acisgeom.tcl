if { [Ng_ACISCommand isACISavailable] == "yes" } {
    .ngmenu.geometry add command -label "ACIS Topology Explorer..." \
        -command { acisdialog; }

    .ngmenu.geometry add command -label "ACIS combine all" \
        -command { Ng_ACISCommand combineall }
    .ngmenu.geometry add command -label "ACIS Create CT" \
        -command { Ng_ACISCommand createct }
}


