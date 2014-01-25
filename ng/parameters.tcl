proc setgranularity { gran } {
#
#    puts "set granularity $gran"
#
    if {$gran == 6} { return }
    set gran [expr $gran - 1]
#
    global options.curvaturesafety
    set surfcurvlist { 1 1.5 2 3 5 }
    set options.curvaturesafety [lindex $surfcurvlist $gran]

    global options.segmentsperedge
    set spelist { 0.3 0.5 1 2 3 }
    set options.segmentsperedge [lindex $spelist $gran]
    
    global stloptions.resthsurfcurvfac
    set surfcurvfaclist { 0.25 0.5 1 1.5 3 }
    set stloptions.resthsurfcurvfac [lindex $surfcurvfaclist $gran]

    global stloptions.resthchartdistfac
    set chartdistfaclist { 0.8 1 1.5 2 5 }
    set stloptions.resthchartdistfac [lindex $chartdistfaclist $gran]

    global stloptions.resthlinelengthfac
    set linelengthfaclist { 0.2 0.35 0.5 1.5 3 }
    set stloptions.resthlinelengthfac [lindex $linelengthfaclist $gran]

    global stloptions.resthcloseedgefac
    set closeedgefaclist { 0.5 1 2 3.5 5 }
    set stloptions.resthcloseedgefac [lindex $closeedgefaclist $gran]

	global stloptions.resthminedgelen
    set minedgelenlist { 0.002 0.02 0.2 1.0 2.0 5.0 10.0 }
    set stloptions.resthminedgelen [lindex $minedgelenlist $gran]
	
    global stloptions.resthedgeanglefac
    set edgeanglefaclist { 0.25 0.5 1 1.5 3 }
    set stloptions.resthedgeanglefac [lindex $edgeanglefaclist $gran]


    global stloptions.resthsurfmeshcurvfac 
    set surfmeshcurvlist { 1 1.5 2 3 5 }
    set stloptions.resthsurfmeshcurvfac [lindex $surfmeshcurvlist $gran]


    global options.grading
    set gradinglist { 0.7 0.5 0.3 0.2 0.1 }
    set options.grading [lindex $gradinglist $gran]
    
}
