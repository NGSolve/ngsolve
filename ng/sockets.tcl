set sockets.serverport 0
set sockets.serverhost "localhost"
set sockets.serverlistbox 0
set sockets.queuelistbox 0
set sockets.currentjoblistbox 0
set sockets.answerlistbox 0
set sockets.myidlabel -1


proc updateserverlist { } {
    global sockets.serverlistbox
    
    set retval [Ng_Socket getserverlist]

    ${sockets.serverlistbox} delete 0 end

    for {set i 0} {$i < [llength $retval]} {incr i 3} {
	${sockets.serverlistbox} insert end \
	    [format "%-16s   %6i   %6i" [lindex $retval $i] [lindex $retval [expr $i+1]] [lindex $retval [expr $i+2]]]
    }
}

proc clientsocketdialog { } {
    set w .clientsock_dlg
    
    if {[winfo exists .clientsock_dlg] == 1} {
	wm withdraw $w
	wm deiconify $w
	focus $w 
    } {
	toplevel $w

	global sockets.serverhost
	global sockets.serverport

	ttk::frame $w.general
	ttk::frame $w.host
	ttk::label $w.host.lab -text "Serverhost: "
	ttk::entry $w.host.name -width 30 -textvariable sockets.serverhost

	pack $w.host.lab $w.host.name -side left
	pack $w.host

	ttk::frame $w.ports
	ttk::label $w.ports.lab1 -text "Serverport: "
	ttk::entry $w.ports.statport -width 6 -textvariable sockets.serverport
	
	pack $w.ports.lab1 $w.ports.statport -side left
	pack $w.ports

	ttk::frame $w.listboxes

	ttk::frame $w.listboxes.choosesocketframe

	tixScrolledListBox $w.listboxes.choosesocketframe.choosesocket -scrollbar auto

	global sockets.serverlistbox

	set sockets.serverlistbox [$w.listboxes.choosesocketframe.choosesocket subwidget listbox]

	${sockets.serverlistbox} configure -width 35
	${sockets.serverlistbox} configure -selectmode browse
	${sockets.serverlistbox} configure -exportselection false

	ttk::button $w.addserver -text "Add ServerSocket" -command {
	    Ng_Socket addserver ${sockets.serverport} ${sockets.serverhost}
	    updateserverlist
	}
	
	pack $w.addserver

	ttk::label $w.linefeed -text "\n"
	pack $w.linefeed
	
	ttk::frame $w.clientidframe
	ttk::label $w.clientidframe.lab -text "Client ID: ";
	global sockets.myidlabel
	ttk::entry $w.clientidframe.val -width 5 -textvariable sockets.myidlabel
	ttk::button $w.clientidframe.but -text "Set" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		Ng_Socket setid $opserver ${sockets.myidlabel}
		updateserverlist
	    }
	}

	pack $w.clientidframe.lab $w.clientidframe.val $w.clientidframe.but -side left
	pack $w.clientidframe


#	label $w.clientidlabel -text "\nClient ID: -1"
#	global sockets.myidlabel
#	set sockets.myidlabel $w.clientidlabel
#	pack $w.clientidlabel


	ttk::label $w.listboxes.choosesocketframe.chooselab -text [format "\n\n%-16s    %6s  %6s                       " Host Socket MyID ]
	pack $w.listboxes.choosesocketframe.chooselab
	pack $w.listboxes.choosesocketframe.choosesocket

	ttk::frame $w.listboxes.choosesocketframe.serverbuttons

	ttk::button $w.listboxes.choosesocketframe.serverbuttons.save -text "Save" -command {
	    Ng_Socket saveserverlist
	}

	global sockets.serverlist
	Ng_Socket loadserverlist
	updateserverlist

	ttk::button $w.listboxes.choosesocketframe.serverbuttons.delete -text "Delete" -command {
	   set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		Ng_Socket deletesocket [lindex $opsel 0]
		updateserverlist
	    } 
	}
	
	pack $w.listboxes.choosesocketframe.serverbuttons.save $w.listboxes.choosesocketframe.serverbuttons.delete -side left
	pack $w.listboxes.choosesocketframe.serverbuttons

	ttk::frame $w.listboxes.statusframe

	ttk::label $w.listboxes.statusframe.statuslabel1 -text "\n\njobqueue"

	tixScrolledListBox $w.listboxes.statusframe.queuestatus -scrollbar auto

	ttk::label $w.listboxes.statusframe.statuslabel2 -text "\ncurrent job"

	tixScrolledListBox $w.listboxes.statusframe.currentjobstatus -scrollbar auto

	ttk::label $w.listboxes.statusframe.statuslabel3 -text "\nanswers"

	tixScrolledListBox $w.listboxes.statusframe.answers -scrollbar auto

	global sockets.queuelistbox
	global sockets.currentjoblistbox
	global sockets.answerlistbox

	set sockets.queuelistbox [$w.listboxes.statusframe.queuestatus subwidget listbox]
	set sockets.currentjoblistbox [$w.listboxes.statusframe.currentjobstatus subwidget listbox]
	set sockets.answerlistbox [$w.listboxes.statusframe.answers subwidget listbox]

	${sockets.queuelistbox} configure -width 50
	${sockets.queuelistbox} configure -height 5
	${sockets.queuelistbox} configure -selectmode browse
	${sockets.queuelistbox} configure -exportselection false
	
	${sockets.currentjoblistbox} configure -width 50
	${sockets.currentjoblistbox} configure -height 1
	${sockets.currentjoblistbox} configure -selectmode browse
	${sockets.currentjoblistbox} configure -exportselection false

	${sockets.answerlistbox} configure -width 50
	${sockets.answerlistbox} configure -height 5
	${sockets.answerlistbox} configure -selectmode browse
	${sockets.answerlistbox} configure -exportselection false

	ttk::button $w.listboxes.statusframe.updatebutton -text "Update" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [Ng_Socket sendqueuestatus $opserver]

		${sockets.queuelistbox} delete 0 end
		
		if {[lindex $retval 0] > 0} {
		    ${sockets.queuelistbox} insert end [format "Blocked for user %i" [lindex $retval 0]]
		} {
		    ${sockets.queuelistbox} insert end "Not blocked"
		}
		
		for {set i 2} {$i < [expr 2*[lindex $retval 1]+2]} {incr i 2} {
		    ${sockets.queuelistbox} insert end [format "client %i, command %s" [lindex $retval $i] [lindex $retval [expr $i+1]]]
		}
		
		${sockets.answerlistbox} delete 0 end
		
		for {set i [expr 2*[lindex $retval 1]+3]} {$i < [llength $retval]} {incr i 2} {
		    ${sockets.answerlistbox} insert end [format "client %i, command %s" [lindex $retval $i] [lindex $retval [expr $i+1]]]
		}

		${sockets.currentjoblistbox} delete 0 end
		set retval [Ng_Socket sendjobstatus $opserver]
		if {[lindex $retval 0] != 0} {
		    ${sockets.currentjoblistbox} insert end [format "client %i, command %s: %s" [lindex $retval 0] [lindex $retval 1] [lrange $retval 2 end]]
		}
		
	    }
	}

	pack $w.listboxes.statusframe.statuslabel1 $w.listboxes.statusframe.queuestatus \
	    $w.listboxes.statusframe.statuslabel2 $w.listboxes.statusframe.currentjobstatus \
	    $w.listboxes.statusframe.statuslabel3 $w.listboxes.statusframe.answers \
	    $w.listboxes.statusframe.updatebutton

	pack $w.listboxes.choosesocketframe $w.listboxes.statusframe -side left
	
	pack $w.listboxes

	ttk::label $w.lab1 -text "\n"
	pack $w.lab1


	ttk::frame $w.buttons1
	ttk::frame $w.buttons2
	
	ttk::button $w.buttons1.getid -text "Get ID" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [Ng_Socket getid $opserver]
		updateserverlist
		set sockets.myidlabel $retval
	    }
	}

	ttk::button $w.buttons1.killjob -text "Kill Cur. Job" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		Ng_Socket killcurrentjob $opserver
	    }
	}

	ttk::button $w.buttons2.sendmesh -text "Send Mesh" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [Ng_Socket sendmesh $opserver]
		set sockets.meshsent 1
	    }
	}

	ttk::button $w.buttons2.sendpde -text "Send PDE" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [NGS_Socket sendpdefile $opserver]
	    }
	}

	ttk::button $w.buttons2.solvepde -text "Solve PDE" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [NGS_Socket solvepde $opserver]
	    }
	}

	ttk::button $w.buttons2.writesol -text "Write Solution" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [NGS_Socket writesolution $opserver]
	    }
	}

	ttk::button $w.buttons2.sendsol -text "Receive Solution" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [NGS_Socket sendsolution $opserver]
	    }
	}

	ttk::button $w.buttons1.blockserver -text "Block Server" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [Ng_Socket blockserver $opserver]
	    }
	}

	ttk::button $w.buttons1.unblockserver -text "UnBlock Server" -command {
	    set opsel [${sockets.serverlistbox} curselection]
	    if {[llength $opsel] > 0} {
		set opserver [lindex $opsel 0]
		set retval [Ng_Socket unblockserver $opserver]
	    }
	}

	
	pack $w.buttons1.getid $w.buttons1.blockserver $w.buttons1.unblockserver $w.buttons1.killjob -side left
	pack $w.buttons2.sendmesh $w.buttons2.sendpde $w.buttons2.solvepde $w.buttons2.writesol $w.buttons2.sendsol -side left

	pack $w.buttons1 $w.buttons2


	wm withdraw $w
	wm geom $w +200+200
	wm deiconify $w
	wm title $w "Client Socket"
	focus .options_dlg

    }
    

}

#.ngmenu.special add command -label "Client Socket" \
    -command { clientsocketdialog }
