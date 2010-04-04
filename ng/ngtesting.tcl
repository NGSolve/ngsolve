
# tests.tcl
proc Ng_TestMeshing { infile outfile logfile} {
  if { ![file exists $infile]} {
    puts $logfile "error: file $infile does not exist"
  } else {
    puts -nonewline $logfile "      loading geometry: ";
    ngtic; Ng_LoadGeometry "$infile"; ngtoc $logfile
    puts -nonewline $logfile "      parsing geometry: ";
    ngtic; Ng_ParseGeometry; ngtoc $logfile
    puts -nonewline $logfile "      generating mesh:  ";
    ngtic; Ng_GenerateMesh; ngtoc $logfile
    puts -nonewline $logfile "      saving mesh:      ";
    ngtic; Ng_SaveMesh "$outfile"; ngtoc $logfile
    flush $logfile
  }
}

# tests.tcl
proc Ngs_TestPDE { infile {nsolves 1} logfile } {
  if { ![file exists $infile]} {
    puts $logfile "error: file $infile does not exist"
  } else {
    puts -nonewline $logfile "       loading PDE file: ";
    ngtic; ngsloadpde "$infile"; ngtoc $logfile
    for {set i 1} {$i<=$nsolves} {incr i 1} {
      puts -nonewline $logfile "      solve PDE level $i: ";
      ngtic; ngssolvepde; ngtoc $logfile
    }
    flush $logfile
  }
}


proc ngtest { {t all} {f ""}} {
  ngset printmsg 0
  if {$f == "" } { 
    set logfile stdout 
  } else {
    set logfile [open $f "w"]
  }
    global options.parthread
    set options.parthread 0
    Ng_SetMeshingParameters 

  if {$t == "all"} {
    ngtest in2d $f
    ngtest geo  $f
    ngtest stl  $f
    ngtest pde  $f
    return
  
  } elseif {$t == "in2d"} {
    puts "\n*** performing in2d file tests ***"
    puts " ** writing results to $f"
    puts "\n ** testing in2d files in examples/ **" 
    set testdir "$::ngdir/../examples"
    set in2dfiles { beam2d hyperbolic piezo2dround rectangle
                    squareincl squareinsquare }
    foreach {tfile} $in2dfiles {
      if {$f != ""} { puts "  * meshing file examples/$tfile.in2d..." }
      puts $logfile "\n  * meshing file examples/$tfile.in2d..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.in2d" "$testdir/$tfile.vol" $logfile
    }  
  
    puts "\n ** testing in2d files in tutorials/ **"
    set testdir "$::ngdir/../share/netgen"
    set in2dfiles { demo2d  newin2d  square  v2in2d }
    foreach {tfile} $in2dfiles {
      if {$f != ""} { puts "  * meshing file tutorials/$tfile.in2d..." }
      puts $logfile "\n  * meshing file tutorials/$tfile.in2d..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.in2d" "$testdir/$tfile.vol" $logfile
    }
    puts "*** in2d tests complete"
  } elseif {$t == "geo"} {
    puts "\n*** performing geo file tests ***"
    puts " ** writing results to $f"
    puts "\n ** testing geo files in examples/ **"
    set testdir "$::ngdir/../examples"
    set geofiles { beam cylsphere rboxcyl thinc 
       boxcyl fichera period saw_3d thinplate 
       coilshield  gamm3d plate shaft tripelpendel 
       cube halfsphere  poly shell twocubes 
       cylinder kaese  quarter3d  skew_prisms }
    foreach {tfile} $geofiles {
      if {$f != ""} { puts "  * meshing file examples/$tfile.geo..."; flush $logfile }
      puts $logfile "\n  * meshing file examples/$tfile.geo..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.geo" "$testdir/$tfile.vol" $logfile
    } 
  
    puts "\n ** testing geo files in tutorials/ **"
    set testdir "$::ngdir/../share/netgen"
    set geofiles { boxcyl          cubemcyl     extrusion  revolution    trafo
                   circle_on_cube  cubemsphere  fichera    sculpture     twobricks
                   cone            cylinder     lshape3d   shaft         twocubes
                   cubeandring     cylsphere    manyholes  sphere        twocyl
                   cubeandspheres  ellipsoid    matrix     sphereincube
                   cube            ellipticcyl  period     torus }
    foreach {tfile} $geofiles {
      if {$f != ""} { puts "  * meshing file $tfile.geo..." }
      puts $logfile "\n  * meshing file $tfile.geo..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.geo" "$testdir/$tfile.vol" $logfile
    } 
  
  } elseif {$t == "stl"} {
    puts "\n*** performing stl file tests ***"  
#    set logfile [open stltest.log "w"]
    puts " ** writing results to $f"
    puts "\n ** testing stl files in examples/ **"
    set testdir "$::ngdir/../examples"
    set stlfiles { crankshaft hinge1 }
    foreach {tfile} $stlfiles {
      if {$f != ""} { puts "  * meshing file examples/$tfile.stl..." }
      puts $logfile "\n  * meshing file examples/$tfile.stl..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.stl" "$testdir/$tfile.vol" $logfile
    }
  
    puts "\n ** testing stl files in tutorials/ **"
    set testdir "$::ngdir/../tutorials"
    set stlfiles { hinge part1 }
    foreach {tfile} $stlfiles {
      if {$f != ""} { puts "  * meshing file tutorials/$tfile.stl..." }
      puts $logfile "\n  * meshing file tutorials/$tfile.stl..."; flush $logfile 
      Ng_TestMeshing "$testdir/$tfile.stl" "$testdir/$tfile.vol" $logfile
    }
    puts "*** stl tests complete"    
  } elseif {$t == "pde"} {
    puts "\n*** preforming pde tests ***"
#    set logfile [open pdetest.log "w"]
    puts " ** writing results to $f"
    
    puts "\n ** testing pde files in ngsolve/pde_tutorial/ **"
    set testdir "$::ngdir/../ngsolve/pde_tutorial"
    set pdefiles { d1_square.pde d2_chip.pde d3_helmholtz.pde d4_cube.pde
                   d5_beam.pde d6_shaft.pde d7_coil.pde d8_coilshield.pde }
     
    foreach {tfile} $pdefiles {
      if {$f != ""} { puts "  * testing ngsolve/pde_tutorial/$tfile..." }
      puts $logfile "\n  * testing ngsolve/pde_tutorial/$tfile..."; flush $logfile 
      Ngs_TestPDE "$testdir/$tfile" 1 $logfile    
    }
    puts "*** pde tests complete"    

  } else {
    puts $logfile "error: unkown test program '$t'"; flush $logfile 
  }
  puts ""
}
Ng_RegisterCmd "ngtest" "netgen" "ngtest" "perform ng standard tests: all in2d geo stl"  




