const char * ngscript[] = {
"\n",\
"\n",\
"if {[catch {package require Tix }]} {\n",\
"    puts \"cannot find package Tix\"\n",\
"}\n",\
"set userlevel 3\n",\
"if { [Ng_GetCommandLineParameter expert]==\"defined\" } {\n",\
"    set userlevel 3\n",\
"}\n",\
"\n",\
"set progname \"NETGEN\"\n",\
"\n",\
"set ngdir \"\"\n",\
"if { [lsearch [array names env] NETGENDIR] != -1 } {\n",\
"    set ngdir $env(NETGENDIR) \n",\
"}\n",\
"if { [string length $ngdir] == 0 } {\n",\
"    set ngdir \".\" \n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"set batchmode [Ng_GetCommandLineParameter batchmode]\n",\
"\n",\
"set solvemode 0\n",\
"if { [Ng_GetCommandLineParameter solve] != \"undefined\" || [Ng_GetCommandLineParameter recent] == \"defined\" } {\n",\
"    set solvemode defined\n",\
"}\n",\
"\n",\
"set shellmode [Ng_GetCommandLineParameter shellmode]\n",\
"\n",\
"if { $shellmode == \"defined\" } {\n",\
"  set batchmode \"defined\"\n",\
"}\n",\
"\n",\
"\n",\
"if { $batchmode != \"defined\" } {\n",\
"	catch {\n",\
"     		wm withdraw .\n",\
"     \n",\
"     		wm title . $progname\n",\
"     		wm geometry . =800x600\n",\
"     		wm minsize . 400 300\n",\
"	}\n",\
"}\n",\
"\n",\
"\n",\
"set drawmode rotate\n",\
"set selectvisual geometry\n",\
"\n",\
"set dirname .\n",\
"set basefilename filename\n",\
"\n",\
"set meshoptions.fineness 3\n",\
"set meshoptions.firststep ag\n",\
"set meshoptions.laststep ov\n",\
"set options.memory 0\n",\
"\n",\
"set options.localh 1\n",\
"set options.delaunay 1\n",\
"set options.checkoverlap 1\n",\
"set options.checkoverlappingboundary 0\n",\
"set options.checkchartboundary 1\n",\
"set options.startinsurface 0\n",\
"set options.blockfill 1\n",\
"set options.debugmode 0\n",\
"set options.dooptimize 1\n",\
"set options.parthread 1\n",\
"set options.elsizeweight 0.2\n",\
"set options.secondorder 0\n",\
"set options.elementorder 1\n",\
"set options.quad 0\n",\
"set options.inverttets 0\n",\
"set options.inverttrigs 0\n",\
"set options.autozrefine 0\n",\
"\n",\
"\n",\
"set options.meshsize 1000\n",\
"set options.minmeshsize 0\n",\
"\n",\
"set options.curvaturesafety 2\n",\
"set options.segmentsperedge 2\n",\
"set options.meshsizefilename \"\"\n",\
"set options.badellimit 175\n",\
"set options.optsteps2d 3\n",\
"set options.optsteps3d 5\n",\
"set options.opterrpow 2\n",\
"\n",\
"set options.grading 0.5\n",\
"set options.printmsg 2\n",\
"\n",\
"set debug.slowchecks 0\n",\
"set debug.debugoutput 0\n",\
"set debug.haltexistingline 0\n",\
"set debug.haltoverlap 0\n",\
"set debug.haltsuccess 0\n",\
"set debug.haltnosuccess 0\n",\
"set debug.haltlargequalclass 0\n",\
"set debug.haltsegment 0\n",\
"set debug.haltnode 0\n",\
"set debug.haltface 0\n",\
"set debug.haltfacenr 0\n",\
"set debug.haltsegmentp1 0\n",\
"set debug.haltsegmentp2 0\n",\
"\n",\
"set geooptions.drawcsg 1\n",\
"set geooptions.detail 0.001\n",\
"set geooptions.accuracy 1e-6\n",\
"set geooptions.facets 20\n",\
"set geooptions.minx -1000\n",\
"set geooptions.miny -1000\n",\
"set geooptions.minz -1000\n",\
"set geooptions.maxx 1000\n",\
"set geooptions.maxy 1000\n",\
"set geooptions.maxz 1000\n",\
"\n",\
"set viewqualityplot 0\n",\
"set memuseplot 0\n",\
"set viewrotatebutton 0\n",\
"set showsensitivehelp 0\n",\
"set showhelpline 0\n",\
"\n",\
"set viewoptions.specpointvlen 0.3\n",\
"set viewoptions.light.amb 0.3\n",\
"set viewoptions.light.diff 0.7\n",\
"set viewoptions.light.spec 1\n",\
"set viewoptions.light.locviewer 0\n",\
"set viewoptions.mat.shininess 50\n",\
"set viewoptions.mat.transp 0.3\n",\
"set viewoptions.colormeshsize 0\n",\
"set viewoptions.whitebackground 1\n",\
"set viewoptions.drawcoordinatecross 1\n",\
"set viewoptions.drawcolorbar 1\n",\
"set viewoptions.drawnetgenlogo 1\n",\
"set viewoptions.stereo 0\n",\
"set viewoptions.shrink 1\n",\
"\n",\
"set viewoptions.drawfilledtrigs 1\n",\
"set viewoptions.drawedges 0\n",\
"set viewoptions.drawbadels 0\n",\
"set viewoptions.centerpoint 0\n",\
"set viewoptions.drawelement 0\n",\
"set viewoptions.drawoutline 1\n",\
"set viewoptions.drawtets 0\n",\
"set viewoptions.drawtetsdomain 0\n",\
"set viewoptions.drawprisms 0\n",\
"set viewoptions.drawpyramids 0\n",\
"set viewoptions.drawhexes 0\n",\
"set viewoptions.drawidentified 0\n",\
"set viewoptions.drawpointnumbers 0\n",\
"set viewoptions.drawedgenumbers 0\n",\
"set viewoptions.drawfacenumbers 0\n",\
"set viewoptions.drawelementnumbers 0\n",\
"set viewoptions.drawdomainsurf 0\n",\
"\n",\
"set viewoptions.drawededges 1\n",\
"set viewoptions.drawedpoints 1\n",\
"set viewoptions.drawedpointnrs 0\n",\
"set viewoptions.drawedtangents 0\n",\
"set viewoptions.drawededgenrs 0\n",\
"set viewoptions.drawmetispartition 0\n",\
"\n",\
"set viewoptions.drawcurveproj 0\n",\
"set viewoptions.drawcurveprojedge 1\n",\
"\n",\
"set viewoptions.clipping.nx 0\n",\
"set viewoptions.clipping.ny 1\n",\
"set viewoptions.clipping.nz 0\n",\
"set viewoptions.clipping.dist 0\n",\
"set viewoptions.clipping.enable 0\n",\
"set viewoptions.clipping.onlydomain 0\n",\
"set viewoptions.clipping.notdomain 0\n",\
"\n",\
"set viewoptions.usecentercoords 0\n",\
"set viewoptions.centerx 0\n",\
"set viewoptions.centery 0\n",\
"set viewoptions.centerz 0\n",\
"\n",\
"set viewoptions.drawspecpoint 0\n",\
"set viewoptions.specpointx 0\n",\
"set viewoptions.specpointy 0\n",\
"set viewoptions.specpointz 0\n",\
"\n",\
"\n",\
"set stloptions.showtrias 0\n",\
"set stloptions.showfilledtrias 1\n",\
"set stloptions.showedges 1\n",\
"set stloptions.showmarktrias 0\n",\
"set stloptions.showactivechart 0\n",\
"set stloptions.yangle 30\n",\
"set stloptions.contyangle 20\n",\
"set stloptions.edgecornerangle 60\n",\
"set stloptions.chartangle 15\n",\
"set stloptions.outerchartangle 70\n",\
"set stloptions.usesearchtree 0\n",\
"set stloptions.chartnumber 1\n",\
"set stloptions.charttrignumber 1\n",\
"set stloptions.chartnumberoffset 0\n",\
"\n",\
"set stloptions.atlasminh 0.1\n",\
"set stloptions.resthsurfcurvfac 2\n",\
"set stloptions.resthsurfcurvenable 0\n",\
"set stloptions.resthatlasfac 2\n",\
"set stloptions.resthatlasenable 1\n",\
"set stloptions.resthchartdistfac 1.2\n",\
"set stloptions.resthchartdistenable 1\n",\
"set stloptions.resthlinelengthfac 0.5\n",\
"set stloptions.resthlinelengthenable 1\n",\
"set stloptions.resthcloseedgefac 1\n",\
"set stloptions.resthcloseedgeenable 1\n",\
"set stloptions.resthedgeanglefac 1\n",\
"set stloptions.resthedgeangleenable 0\n",\
"set stloptions.resthsurfmeshcurvfac 1\n",\
"set stloptions.resthsurfmeshcurvenable 0\n",\
"set stloptions.recalchopt 1\n",\
"\n",\
"set stldoctor.drawmeshededges 1\n",\
"set stldoctor.geom_tol_fact 0.000001\n",\
"set stldoctor.useexternaledges 0\n",\
"set stldoctor.showfaces 0\n",\
"set stldoctor.conecheck 1\n",\
"set stldoctor.spiralcheck 1\n",\
"set stldoctor.selecttrig 0\n",\
"set stldoctor.selectmode 1\n",\
"set stldoctor.longlinefact 0\n",\
"set stldoctor.showexcluded 1\n",\
"set stldoctor.edgeselectmode 0\n",\
"set stldoctor.nodeofseltrig 1\n",\
"set stldoctor.showtouchedtrigchart 0\n",\
"set stldoctor.showedgecornerpoints 0\n",\
"set stldoctor.showmarkedtrigs 1\n",\
"set stldoctor.dirtytrigfact 0.01\n",\
"set stldoctor.smoothangle 90\n",\
"set stldoctor.selectwithmouse 1\n",\
"set stldoctor.showvicinity 0\n",\
"set stldoctor.vicinity 50\n",\
"set stldoctor.smoothnormalsweight 0.2\n",\
"\n",\
"set occoptions.showvolumenr 0\n",\
"set occoptions.showsurfaces 1\n",\
"set occoptions.showedges 1\n",\
"set occoptions.showsolidnr 0\n",\
"set occoptions.showsolidnr2 0\n",\
"set occoptions.visproblemfaces 0\n",\
"set occoptions.zoomtohighlightedentity 0\n",\
"set occoptions.deflection 1\n",\
"set occoptions.tolerance 1e-3\n",\
"set occoptions.fixsmalledges 1\n",\
"set occoptions.fixspotstripfaces 1\n",\
"set occoptions.sewfaces 1\n",\
"set occoptions.makesolids 1\n",\
"set occoptions.splitpartitions 0\n",\
"\n",\
"set meshdoctor.active 0\n",\
"set meshdoctor.markedgedist 1\n",\
"\n",\
"\n",\
"set status_np 0\n",\
"set status_ne 0\n",\
"set status_nse 0\n",\
"set status_working \" \"\n",\
"set status_task \" \"\n",\
"set status_percent 0\n",\
"set status_filename 0\n",\
"set status_tetqualclasses \"10 20 30 40 10 20 30 40 10 20 30 40 10 20 30 40 10 20 30 40\"\n",\
"\n",\
"set exportfiletype PERMAS\n",\
"\n",\
"set preproc.facenr 0\n",\
"set preproc.selectmode query\n",\
"set preproc.numtrig 0\n",\
"\n",\
"set mem_moveable 0\n",\
"\n",\
"\n",\
"set multithread_pause 0\n",\
"set multithread_testmode 0\n",\
"set multithread_redraw 0\n",\
"set multithread_drawing 0\n",\
"set multithread_terminate 0\n",\
"set multithread_running 0\n",\
"\n",\
"set level 0\n",\
"\n",\
"\n",\
"set tablesforoutput {}\n",\
"\n",\
"\n",\
"\n",\
"set optlist {\n",\
"    options.localh \n",\
"    options.delaunay \n",\
"    options.checkoverlap \n",\
"    options.startinsurface \n",\
"    options.blockfill \n",\
"    options.dooptimize \n",\
"    options.elsizeweight \n",\
"    options.meshsize \n",\
"    options.minmeshsize \n",\
"    options.curvaturesafety \n",\
"    options.optsteps2d \n",\
"    options.optsteps3d \n",\
"    options.secondorder\n",\
"}\n",\
"\n",\
"\n",\
"set visoptions.usetexture 0\n",\
"set visoptions.invcolor 0\n",\
"set visoptions.imaginary 0\n",\
"set visoptions.lineartexture 1\n",\
"set visoptions.numtexturecols 16\n",\
"set visoptions.showclipsolution 1\n",\
"set visoptions.showsurfacesolution 0\n",\
"set visoptions.drawfieldlines 0\n",\
"set visoptions.drawpointcurves 1\n",\
"set visoptions.numfieldlines 100\n",\
"set visoptions.fieldlinesrandomstart 0\n",\
"set visoptions.fieldlinesstartarea box\n",\
"set visoptions.fieldlinesstartareap1x 1\n",\
"set visoptions.fieldlinesstartareap1y 1\n",\
"set visoptions.fieldlinesstartareap1z 1\n",\
"set visoptions.fieldlinesstartareap2x 0\n",\
"set visoptions.fieldlinesstartareap2y 0\n",\
"set visoptions.fieldlinesstartareap2z 0\n",\
"set visoptions.fieldlinesstartface -1\n",\
"set visoptions.fieldlinesfilename none\n",\
"set visoptions.fieldlinestolerance 0.0005\n",\
"set visoptions.fieldlinesrktype crungekutta\n",\
"set visoptions.fieldlineslength 0.5\n",\
"set visoptions.fieldlinesmaxpoints 500\n",\
"set visoptions.fieldlinesthickness 0.0015\n",\
"set visoptions.fieldlinesvecfunction none\n",\
"set visoptions.fieldlinesphase 0\n",\
"set visoptions.fieldlinesonlyonephase 1\n",\
"\n",\
"\n",\
"set visoptions.lineplotfile empty\n",\
"set visoptions.lineplotsource file\n",\
"set visoptions.lineplotusingx 0\n",\
"set visoptions.lineplotusingy 1\n",\
"set visoptions.lineplotautoscale 1\n",\
"set visoptions.lineplotxmin 0\n",\
"set visoptions.lineplotxmax 1\n",\
"set visoptions.lineplotymin 0\n",\
"set visoptions.lineplotymax 1\n",\
"set visoptions.lineplotcurrentnum -1\n",\
"set visoptions.lineplotinfos \"\"\n",\
"set visoptions.lineplotselected none\n",\
"set visoptions.lineplotselector \"\"\n",\
"set visoptions.lineplotcolor red\n",\
"set visoptions.lineplotsizex 500\n",\
"set visoptions.lineplotsizey 400\n",\
"set visoptions.lineplotselectedeval 0\n",\
"set visoptions.lineplotdatadescr \"column1 column2 column3\"\n",\
"set visoptions.lineplotxcoordselector \"\"\n",\
"set visoptions.lineplotycoordselector \"\"\n",\
"set visoptions.evaluatefilenames none\n",\
"set visoptions.evaluatefiledescriptions none\n",\
"\n",\
"\n",\
"set visoptions.clipsolution none\n",\
"set visoptions.scalfunction none\n",\
"set visoptions.vecfunction none\n",\
"set visoptions.evaluate abs\n",\
"set visoptions.gridsize 20\n",\
"set visoptions.xoffset 0\n",\
"set visoptions.yoffset 0\n",\
"set visoptions.autoscale 1\n",\
"set visoptions.lineartexture 1\n",\
"set visoptions.redrawperiodic 0\n",\
"set visoptions.logscale 0\n",\
"set visoptions.mminval 0\n",\
"set visoptions.mmaxval 1\n",\
"set visoptions.isolines 0\n",\
"set visoptions.isosurf 0\n",\
"set visoptions.subdivisions 1\n",\
"set visoptions.numiso 10\n",\
"set visoptions.autoredraw 0\n",\
"set visoptions.autoredrawtime 2\n",\
"set visoptions.simulationtime 0\n",\
"set visoptions.multidimcomponent 0\n",\
"\n",\
"set visoptions.deformation 0\n",\
"set visoptions.scaledeform1 1\n",\
"set visoptions.scaledeform2 1\n",\
"\n",\
"set parallel_netgen 0\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set optfilename ng.opt\n",\
"global env\n",\
"if { [llength [array names env NG_OPT]] == 1 } {\n",\
"    if { [string length $env(NG_OPT)] > 0 } {\n",\
"	set optfilename $env(NG_OPT) \n",\
"    }\n",\
"}\n",\
"\n",\
"if { [file exists $optfilename] == 1 } {\n",\
"    set datei [open $optfilename r]\n",\
"    while { [gets $datei line] >= 0 } {\n",\
"	set [lindex $line 0] [lindex $line 1]\n",\
"    }\n",\
"    close $datei\n",\
"} {\n",\
"    puts \"optfile $optfilename does not exist - using default values\"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc saveoptions { } {\n",\
"    uplevel 1  {\n",\
"	set file ng.opt \n",\
"	\n",\
"	if {$file != \"\"} {\n",\
"	    set datei [open $file w]\n",\
"	    \n",\
"	    puts $datei \"meshoptions.fineness  ${meshoptions.fineness}\"\n",\
"	    puts $datei \"meshoptions.firststep ${meshoptions.firststep}\"\n",\
"	    puts $datei \"meshoptions.laststep  ${meshoptions.laststep}\" \n",\
"	    puts $datei \"options.localh  ${options.localh}\"\n",\
"	    puts $datei \"options.delaunay  ${options.delaunay}\"\n",\
"	    puts $datei \"options.checkoverlap  ${options.checkoverlap}\"\n",\
"	    puts $datei \"options.checkchartboundary  ${options.checkchartboundary}\"\n",\
"	    puts $datei \"options.startinsurface  ${options.startinsurface}\" \n",\
"	    puts $datei \"options.blockfill  ${options.blockfill}\" \n",\
"	    puts $datei \"options.debugmode  ${options.debugmode}\" \n",\
"	    puts $datei \"options.dooptimize ${options.dooptimize}\" \n",\
"	    puts $datei \"options.parthread  ${options.parthread}\"  \n",\
"	    puts $datei \"options.elsizeweight  ${options.elsizeweight}\" \n",\
"	    puts $datei \"options.secondorder  ${options.secondorder}\" \n",\
"	    puts $datei \"options.elementorder  ${options.elementorder}\" \n",\
"	    puts $datei \"options.memory  ${options.memory}\" \n",\
"	    puts $datei \"options.quad  ${options.quad}\" \n",\
"	    puts $datei \"options.inverttets  ${options.inverttets}\" \n",\
"	    puts $datei \"options.inverttrigs  ${options.inverttrigs}\" \n",\
"	    puts $datei \"options.autozrefine ${options.autozrefine}\" \n",\
"	    puts $datei \"options.meshsize  ${options.meshsize}\" \n",\
"	    puts $datei \"options.minmeshsize  ${options.minmeshsize}\" \n",\
"	    puts $datei \"options.curvaturesafety  ${options.curvaturesafety}\" \n",\
"	    puts $datei \"options.segmentsperedge  ${options.segmentsperedge}\" \n",\
"	    puts $datei \"options.meshsizefilename  ${options.meshsizefilename}\" \n",\
"	    puts $datei \"options.badellimit  ${options.badellimit}\" \n",\
"	    puts $datei \"options.optsteps2d  ${options.optsteps2d}\" \n",\
"	    puts $datei \"options.optsteps3d  ${options.optsteps3d}\" \n",\
"	    puts $datei \"options.opterrpow  ${options.opterrpow}\" \n",\
"	    puts $datei \"options.grading  ${options.grading}\" \n",\
"	    puts $datei \"options.printmsg  ${options.printmsg}\" \n",\
"	    puts $datei \"geooptions.drawcsg  ${geooptions.drawcsg}\" \n",\
"	    puts $datei \"geooptions.detail  ${geooptions.detail}\" \n",\
"	    puts $datei \"geooptions.accuracy  ${geooptions.accuracy}\" \n",\
"	    puts $datei \"geooptions.facets  ${geooptions.facets}\" \n",\
"	    puts $datei \"geooptions.minx  ${geooptions.minx}\" \n",\
"	    puts $datei \"geooptions.miny  ${geooptions.miny}\" \n",\
"	    puts $datei \"geooptions.minz  ${geooptions.minz}\" \n",\
"	    puts $datei \"geooptions.maxx  ${geooptions.maxx}\" \n",\
"	    puts $datei \"geooptions.maxy  ${geooptions.maxy}\" \n",\
"	    puts $datei \"geooptions.maxz  ${geooptions.maxz}\" \n",\
"	    puts $datei \"viewoptions.specpointvlen  ${viewoptions.specpointvlen}\" \n",\
"	    puts $datei \"viewoptions.light.amb  ${viewoptions.light.amb}\" \n",\
"	    puts $datei \"viewoptions.light.diff ${viewoptions.light.diff}\"\n",\
"	    puts $datei \"viewoptions.light.spec ${viewoptions.light.spec}\"\n",\
"	    puts $datei \"viewoptions.light.locviewer ${viewoptions.light.locviewer}\"\n",\
"	    puts $datei \"viewoptions.mat.shininess  ${viewoptions.mat.shininess}\" \n",\
"	    puts $datei \"viewoptions.mat.transp  ${viewoptions.mat.transp}\" \n",\
"	    puts $datei \"viewoptions.colormeshsize ${viewoptions.colormeshsize}\"\n",\
"	    puts $datei \"viewoptions.whitebackground  ${viewoptions.whitebackground}\" \n",\
"	    puts $datei \"viewoptions.drawcolorbar  ${viewoptions.drawcolorbar}\" \n",\
"	    puts $datei \"viewoptions.drawcoordinatecross  ${viewoptions.drawcoordinatecross}\" \n",\
"	    puts $datei \"viewoptions.drawnetgenlogo  ${viewoptions.drawnetgenlogo}\" \n",\
"	    puts $datei \"viewoptions.stereo  ${viewoptions.stereo}\" \n",\
"	    puts $datei \"viewoptions.drawfilledtrigs  ${viewoptions.drawfilledtrigs}\" \n",\
"	    puts $datei \"viewoptions.drawedges  ${viewoptions.drawedges}\" \n",\
"	    puts $datei \"viewoptions.drawbadels  ${viewoptions.drawbadels}\" \n",\
"	    puts $datei \"viewoptions.centerpoint  ${viewoptions.centerpoint}\" \n",\
"	    puts $datei \"viewoptions.drawelement  ${viewoptions.drawelement}\" \n",\
"	    puts $datei \"viewoptions.drawoutline  ${viewoptions.drawoutline}\" \n",\
"	    puts $datei \"viewoptions.drawtets  ${viewoptions.drawtets}\"\n",\
"	    puts $datei \"viewoptions.drawprisms  ${viewoptions.drawprisms}\"\n",\
"	    puts $datei \"viewoptions.drawpyramids  ${viewoptions.drawpyramids}\" \n",\
"	    puts $datei \"viewoptions.drawhexes  ${viewoptions.drawhexes}\" \n",\
"	    puts $datei \"viewoptions.drawidentified  ${viewoptions.drawidentified}\" \n",\
"	    puts $datei \"viewoptions.drawpointnumbers  ${viewoptions.drawpointnumbers}\" \n",\
"	    \n",\
"	    puts $datei \"viewoptions.drawededges  ${viewoptions.drawededges}\" \n",\
"	    puts $datei \"viewoptions.drawedpoints  ${viewoptions.drawedpoints}\" \n",\
"	    puts $datei \"viewoptions.drawedpointnrs  ${viewoptions.drawedpointnrs}\" \n",\
"	    puts $datei \"viewoptions.drawedtangents  ${viewoptions.drawedtangents}\" \n",\
"	    puts $datei \"viewoptions.shrink  ${viewoptions.shrink}\" \n",\
"	    \n",\
"	    puts $datei \"stloptions.showtrias  ${stloptions.showtrias}\" \n",\
"	    puts $datei \"stloptions.showfilledtrias  ${stloptions.showfilledtrias}\" \n",\
"	    puts $datei \"stloptions.showedges  ${stloptions.showedges}\" \n",\
"	    puts $datei \"stloptions.showmarktrias  ${stloptions.showmarktrias}\" \n",\
"	    puts $datei \"stloptions.showactivechart  ${stloptions.showactivechart}\" \n",\
"	    puts $datei \"stloptions.yangle  ${stloptions.yangle}\" \n",\
"	    puts $datei \"stloptions.contyangle  ${stloptions.contyangle}\" \n",\
"	    puts $datei \"stloptions.edgecornerangle  ${stloptions.edgecornerangle}\" \n",\
"	    puts $datei \"stloptions.chartangle  ${stloptions.chartangle}\" \n",\
"	    puts $datei \"stloptions.outerchartangle  ${stloptions.outerchartangle}\" \n",\
"	    puts $datei \"stloptions.usesearchtree  ${stloptions.usesearchtree}\" \n",\
"	    puts $datei \"stloptions.chartnumber  ${stloptions.chartnumber}\" \n",\
"	    puts $datei \"stloptions.charttrignumber  ${stloptions.charttrignumber}\" \n",\
"	    puts $datei \"stloptions.chartnumberoffset  ${stloptions.chartnumberoffset}\" \n",\
"	    puts $datei \"stloptions.atlasminh  ${stloptions.atlasminh}\" \n",\
"	    puts $datei \"stloptions.resthsurfcurvfac  ${stloptions.resthsurfcurvfac}\" \n",\
"	    puts $datei \"stloptions.resthsurfcurvenable  ${stloptions.resthsurfcurvenable}\" \n",\
"	    puts $datei \"stloptions.resthatlasfac  ${stloptions.resthatlasfac}\" \n",\
"	    puts $datei \"stloptions.resthatlasenable  ${stloptions.resthatlasenable}\" \n",\
"	    puts $datei \"stloptions.resthchartdistfac  ${stloptions.resthchartdistfac}\" \n",\
"	    puts $datei \"stloptions.resthchartdistenable  ${stloptions.resthchartdistenable}\" \n",\
"	    puts $datei \"stloptions.resthlinelengthfac  ${stloptions.resthlinelengthfac}\" \n",\
"	    puts $datei \"stloptions.resthlinelengthenable  ${stloptions.resthlinelengthenable}\" \n",\
"	    puts $datei \"stloptions.resthcloseedgefac  ${stloptions.resthcloseedgefac}\" \n",\
"	    puts $datei \"stloptions.resthcloseedgeenable  ${stloptions.resthcloseedgeenable}\" \n",\
"	    puts $datei \"stloptions.resthedgeanglefac  ${stloptions.resthedgeanglefac}\" \n",\
"	    puts $datei \"stloptions.resthedgeangleenable  ${stloptions.resthedgeangleenable}\" \n",\
"	    puts $datei \"stloptions.resthsurfmeshcurvfac  ${stloptions.resthsurfmeshcurvfac}\" \n",\
"	    puts $datei \"stloptions.resthsurfmeshcurvenable  ${stloptions.resthsurfmeshcurvenable}\" \n",\
"	    puts $datei \"stloptions.recalchopt  ${stloptions.recalchopt}\" \n",\
"	    \n",\
"	    puts $datei \"visoptions.subdivisions ${visoptions.subdivisions}\"\n",\
"\n",\
"\n",\
"	    	    	    if { [info exists trafooptions.solver] == 1 } {\n",\
"		puts $datei \"trafooptions.solver ${trafooptions.solver}\" \n",\
"		puts $datei \"trafooptions.levels ${trafooptions.levels}\" \n",\
"		puts $datei \"trafooptions.linits ${trafooptions.linits}\" \n",\
"		puts $datei \"trafooptions.nonlinits ${trafooptions.nonlinits}\" \n",\
"		puts $datei \"trafooptions.stabcurrent ${trafooptions.stabcurrent}\" \n",\
"		puts $datei \"trafooptions.checkcond ${trafooptions.checkcond}\" \n",\
"		puts $datei \"trafooptions.maxdirect ${trafooptions.maxdirect}\" \n",\
"		puts $datei \"trafooptions.secondorder ${trafooptions.secondorder}\" \n",\
"		puts $datei \"trafooptions.homogenizedcore ${trafooptions.homogenizedcore}\" \n",\
"		puts $datei \"trafooptions.ordercore ${trafooptions.ordercore}\" \n",\
"		puts $datei \"trafooptions.simplecurrents ${trafooptions.simplecurrents}\" \n",\
"		puts $datei \"trafooptions.assemblecomplexmatrix ${trafooptions.assemblecomplexmatrix}\" \n",\
"\n",\
"		puts $datei \"trafooptions.meshcasing  ${trafooptions.meshcasing}\" \n",\
"		puts $datei \"trafooptions.meshcore    ${trafooptions.meshcore}\" \n",\
"		puts $datei \"trafooptions.meshclumps  ${trafooptions.meshclumps}\" \n",\
"		puts $datei \"trafooptions.meshshields ${trafooptions.meshshields}\" \n",\
"		puts $datei \"trafooptions.meshcoils   ${trafooptions.meshcoils}\" \n",\
"		puts $datei \"trafooptions.bcmdirectory  ${trafooptions.bcmdirectory}\" \n",\
"		puts $datei \"trafooptions.lossdensityfile  ${trafooptions.lossdensityfile}\" \n",\
"	    }\n",\
"\n",\
"	    if { [info exists smalltrafomodell.tankheight] == 1 } {\n",\
"		puts $datei \"smalltrafomodell.tankheight ${smalltrafomodell.tankheight}\"\n",\
"		puts $datei \"smalltrafomodell.tankwidth ${smalltrafomodell.tankwidth}\"\n",\
"		puts $datei \"smalltrafomodell.tanklength ${smalltrafomodell.tanklength}\"\n",\
"		puts $datei \"smalltrafomodell.corewidth ${smalltrafomodell.corewidth}\"\n",\
"		puts $datei \"smalltrafomodell.windowheight ${smalltrafomodell.windowheight}\"\n",\
"		puts $datei \"smalltrafomodell.limbdistance ${smalltrafomodell.limbdistance}\"\n",\
"		puts $datei \"smalltrafomodell.xposcore ${smalltrafomodell.xposcore}\"\n",\
"		puts $datei \"smalltrafomodell.yposcore ${smalltrafomodell.yposcore}\"\n",\
"		puts $datei \"smalltrafomodell.zposcore ${smalltrafomodell.zposcore}\"\n",\
"		puts $datei \"smalltrafomodell.leakagefluxguidethickness ${smalltrafomodell.leakagefluxguidethickness}\"\n",\
"		puts $datei \"smalltrafomodell.leakagefluxguidewidth ${smalltrafomodell.leakagefluxguidewidth}\"\n",\
"		puts $datei \"smalltrafomodell.leakagefluxguidezposition ${smalltrafomodell.leakagefluxguidezposition}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.1 ${smalltrafomodell.limbcoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.1 ${smalltrafomodell.ricoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.1 ${smalltrafomodell.rocoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.1 ${smalltrafomodell.zposcoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.1 ${smalltrafomodell.heightcoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.1 ${smalltrafomodell.currentcoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.1 ${smalltrafomodell.nturnscoil.1}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.2 ${smalltrafomodell.limbcoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.2 ${smalltrafomodell.ricoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.2 ${smalltrafomodell.rocoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.2 ${smalltrafomodell.zposcoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.2 ${smalltrafomodell.heightcoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.2 ${smalltrafomodell.currentcoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.2 ${smalltrafomodell.nturnscoil.2}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.3 ${smalltrafomodell.limbcoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.3 ${smalltrafomodell.ricoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.3 ${smalltrafomodell.rocoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.3 ${smalltrafomodell.zposcoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.3 ${smalltrafomodell.heightcoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.3 ${smalltrafomodell.currentcoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.3 ${smalltrafomodell.nturnscoil.3}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.4 ${smalltrafomodell.limbcoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.4 ${smalltrafomodell.ricoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.4 ${smalltrafomodell.rocoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.4 ${smalltrafomodell.zposcoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.4 ${smalltrafomodell.heightcoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.4 ${smalltrafomodell.currentcoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.4 ${smalltrafomodell.nturnscoil.4}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.5 ${smalltrafomodell.limbcoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.5 ${smalltrafomodell.ricoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.5 ${smalltrafomodell.rocoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.5 ${smalltrafomodell.zposcoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.5 ${smalltrafomodell.heightcoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.5 ${smalltrafomodell.currentcoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.5 ${smalltrafomodell.nturnscoil.5}\"\n",\
"		puts $datei \"smalltrafomodell.limbcoil.6 ${smalltrafomodell.limbcoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.ricoil.6 ${smalltrafomodell.ricoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.rocoil.6 ${smalltrafomodell.rocoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.zposcoil.6 ${smalltrafomodell.zposcoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.heightcoil.6 ${smalltrafomodell.heightcoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.currentcoil.6 ${smalltrafomodell.currentcoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.nturnscoil.6 ${smalltrafomodell.nturnscoil.6}\"\n",\
"		puts $datei \"smalltrafomodell.limbtest.1 ${smalltrafomodell.limbtest.1}\"\n",\
"		puts $datei \"smalltrafomodell.heighttest.1 ${smalltrafomodell.heighttest.1}\"\n",\
"		puts $datei \"smalltrafomodell.widthtest.1 ${smalltrafomodell.widthtest.1}\"\n",\
"		puts $datei \"smalltrafomodell.rtest.1 ${smalltrafomodell.rtest.1}\"\n",\
"		puts $datei \"smalltrafomodell.zpostest.1 ${smalltrafomodell.zpostest.1}\"\n",\
"		puts $datei \"smalltrafomodell.edgeradiustest.1 ${smalltrafomodell.edgeradiustest.1}\"\n",\
"		puts $datei \"smalltrafomodell.finetest.1 ${smalltrafomodell.finetest.1}\"\n",\
"		puts $datei \"smalltrafomodell.conductivetest.1 ${smalltrafomodell.conductivetest.1}\"\n",\
"		puts $datei \"smalltrafomodell.limbtest.2 ${smalltrafomodell.limbtest.2}\"\n",\
"		puts $datei \"smalltrafomodell.heighttest.2 ${smalltrafomodell.heighttest.2}\"\n",\
"		puts $datei \"smalltrafomodell.widthtest.2 ${smalltrafomodell.widthtest.2}\"\n",\
"		puts $datei \"smalltrafomodell.rtest.2 ${smalltrafomodell.rtest.2}\"\n",\
"		puts $datei \"smalltrafomodell.zpostest.2 ${smalltrafomodell.zpostest.2}\"\n",\
"		puts $datei \"smalltrafomodell.edgeradiustest.2 ${smalltrafomodell.edgeradiustest.2}\"\n",\
"		puts $datei \"smalltrafomodell.finetest.2 ${smalltrafomodell.finetest.2}\"\n",\
"		puts $datei \"smalltrafomodell.conductivetest.2 ${smalltrafomodell.conductivetest.2}\"\n",\
"		puts $datei \"smalltrafomodell.limbtest.3 ${smalltrafomodell.limbtest.3}\"\n",\
"		puts $datei \"smalltrafomodell.heighttest.3 ${smalltrafomodell.heighttest.3}\"\n",\
"		puts $datei \"smalltrafomodell.widthtest.3 ${smalltrafomodell.widthtest.3}\"\n",\
"		puts $datei \"smalltrafomodell.rtest.3 ${smalltrafomodell.rtest.3}\"\n",\
"		puts $datei \"smalltrafomodell.zpostest.3 ${smalltrafomodell.zpostest.3}\"\n",\
"		puts $datei \"smalltrafomodell.edgeradiustest.3 ${smalltrafomodell.edgeradiustest.3}\"\n",\
"		puts $datei \"smalltrafomodell.finetest.3 ${smalltrafomodell.finetest.3}\"\n",\
"		puts $datei \"smalltrafomodell.conductivetest.3 ${smalltrafomodell.conductivetest.3}\"\n",\
"		puts $datei \"smalltrafomodell.limbtest.4 ${smalltrafomodell.limbtest.4}\"\n",\
"		puts $datei \"smalltrafomodell.heighttest.4 ${smalltrafomodell.heighttest.4}\"\n",\
"		puts $datei \"smalltrafomodell.widthtest.4 ${smalltrafomodell.widthtest.4}\"\n",\
"		puts $datei \"smalltrafomodell.rtest.4 ${smalltrafomodell.rtest.4}\"\n",\
"		puts $datei \"smalltrafomodell.zpostest.4 ${smalltrafomodell.zpostest.4}\"\n",\
"		puts $datei \"smalltrafomodell.edgeradiustest.4 ${smalltrafomodell.edgeradiustest.4}\"\n",\
"		puts $datei \"smalltrafomodell.finetest.4 ${smalltrafomodell.finetest.4}\"\n",\
"		puts $datei \"smalltrafomodell.conductivetest.4 ${smalltrafomodell.conductivetest.4}\"\n",\
"		puts $datei \"smalltrafomodell.nperitest ${smalltrafomodell.nperitest}\"\n",\
"		puts $datei \"smalltrafomodell.filename ${smalltrafomodell.filename}\"\n",\
"		puts $datei \"smalltrafomodell.murlfguide ${smalltrafomodell.murlfguide}\"\n",\
"		puts $datei \"smalltrafomodell.murtestwire ${smalltrafomodell.murtestwire}\"\n",\
"		puts $datei \"smalltrafomodell.murcore ${smalltrafomodell.murcore}\"\n",\
"		puts $datei \"smalltrafomodell.kappalfguide ${smalltrafomodell.kappalfguide}\"\n",\
"		puts $datei \"smalltrafomodell.kappatestwire ${smalltrafomodell.kappatestwire}\"\n",\
"		puts $datei \"smalltrafomodell.kappacore ${smalltrafomodell.kappacore}\"\n",\
"	    }\n",\
"	    \n",\
"	    \n",\
"	    close $datei\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc saveinifile { } {\n",\
"    uplevel 1  {\n",\
"	set datei [open ng.ini w]\n",\
"	for { set i [.ngmenu.file.recent index last] } { $i >= 1 } { incr i -1 } {\n",\
"	    puts $datei \"recentfile \\\"[.ngmenu.file.recent entrycget $i -label]\\\"\"\n",\
"	}\n",\
"	\n",\
"	close $datei\n",\
"    }    \n",\
"\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc savemeshinifile { } {\n",\
"    uplevel 1  {\n",\
"	set datei [open ngmesh.ini w]\n",\
"	for { set i [.ngmenu.file.recentmesh index last] } { $i >= 1 } { incr i -1 } {\n",\
"	    puts $datei \"recentfile \\\"[.ngmenu.file.recentmesh entrycget $i -label]\\\"\"\n",\
"	}\n",\
"	\n",\
"	close $datei\n",\
"    }    \n",\
"\n",\
"\n",\
"}\n",\
"\n",\
"proc loadinifile { } {\n",\
"    if { [file exists ng.ini] == 1 } {\n",\
"	set datei [open ng.ini r]\n",\
"	while { [gets $datei line] >= 0 } {\n",\
"	    if {[lindex $line 0] == \"recentfile\"} {\n",\
"		set filename [lindex $line 1]\n",\
"		AddRecentFile $filename\n",\
"	    }\n",\
"	}\n",\
"	close $datei\n",\
"    }\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"proc loadmeshinifile { } {\n",\
"    if { [file exists ngmesh.ini] == 1 } {\n",\
"	set datei [open ngmesh.ini r]\n",\
"	while { [gets $datei line] >= 0 } {\n",\
"	    if {[lindex $line 0] == \"recentfile\"} {\n",\
"		set filename [lindex $line 1]\n",\
"		AddRecentMeshFile $filename\n",\
"	    }\n",\
"	}\n",\
"	close $datei\n",\
"    }\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"proc setgranularity { gran } {\n",\
"    if {$gran == 6} { return }\n",\
"    set gran [expr $gran - 1]\n",\
"    global options.curvaturesafety\n",\
"    set surfcurvlist { 1 1.5 2 3 5 }\n",\
"    set options.curvaturesafety [lindex $surfcurvlist $gran]\n",\
"\n",\
"    global options.segmentsperedge\n",\
"    set spelist { 0.3 0.5 1 2 3 }\n",\
"    set options.segmentsperedge [lindex $spelist $gran]\n",\
"    \n",\
"    global stloptions.resthsurfcurvfac\n",\
"    set surfcurvfaclist { 0.25 0.5 1 1.5 3 }\n",\
"    set stloptions.resthsurfcurvfac [lindex $surfcurvfaclist $gran]\n",\
"\n",\
"    global stloptions.resthchartdistfac\n",\
"    set chartdistfaclist { 0.8 1 1.5 2 5 }\n",\
"    set stloptions.resthchartdistfac [lindex $chartdistfaclist $gran]\n",\
"\n",\
"    global stloptions.resthlinelengthfac\n",\
"    set linelengthfaclist { 0.2 0.35 0.5 1.5 3 }\n",\
"    set stloptions.resthlinelengthfac [lindex $linelengthfaclist $gran]\n",\
"\n",\
"    global stloptions.resthcloseedgefac\n",\
"    set closeedgefaclist { 0.5 1 2 3.5 5 }\n",\
"    set stloptions.resthcloseedgefac [lindex $closeedgefaclist $gran]\n",\
"\n",\
"    global stloptions.resthedgeanglefac\n",\
"    set edgeanglefaclist { 0.25 0.5 1 1.5 3 }\n",\
"    set stloptions.resthedgeanglefac [lindex $edgeanglefaclist $gran]\n",\
"\n",\
"\n",\
"    global stloptions.resthsurfmeshcurvfac \n",\
"    set surfmeshcurvlist { 1 1.5 2 3 5 }\n",\
"    set stloptions.resthsurfmeshcurvfac [lindex $surfmeshcurvlist $gran]\n",\
"\n",\
"\n",\
"    global options.grading\n",\
"    set gradinglist { 0.7 0.5 0.3 0.2 0.1 }\n",\
"    set options.grading [lindex $gradinglist $gran]\n",\
"    \n",\
"}\n",\
"\n",\
"\n",\
"if { $batchmode != \"defined\" } {\n",\
"    \n",\
"\n",\
"menu .ngmenu -tearoff 0  -relief raised -bd 2\n",\
". configure -menu .ngmenu\n",\
"\n",\
".ngmenu add cascade -label \"File\" -menu .ngmenu.file -underline 0\n",\
".ngmenu add cascade -label \"Geometry\" -menu .ngmenu.geometry -underline 0\n",\
".ngmenu add cascade -label \"Mesh\" -menu .ngmenu.mesh -underline 0\n",\
".ngmenu add cascade -label \"View\" -menu .ngmenu.view -underline 0\n",\
".ngmenu add cascade -label \"Refinement\" -menu .ngmenu.meshsize -underline 5\n",\
"\n",\
"if { $userlevel == 3} {\n",\
"    .ngmenu add cascade -label \"Special\" -menu .ngmenu.special -underline 3\n",\
"}\n",\
"\n",\
".ngmenu add cascade -label \"Help\" -menu .ngmenu.help -underline 0\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.file\n",\
"\n",\
".ngmenu.file add command -label \"Load Geometry...\" -accelerator \"<l><g>\" \\\n",\
"    -command { \n",\
"	set types {\n",\
"	    {\"All Geometry types\"   { .stl .stlb .step .stp .geo .in2d .igs .iges .brep .in2dnew .sat} }\n",\
"	    {\"IGES Geometry\"	{.igs .iges} }\n",\
"	    {\"BREP OpenCascade Geometry\"    {.brep} }\n",\
"	    {\"STL Geometry\"        {.stl} }\n",\
"	    {\"Binary STL Geometry\"    {.stlb} }\n",\
"	    {\"STEP Geometry\"    {.step .stp} }\n",\
"	    {\"Geometry file\"       {.geo} }\n",\
"	    {\"2D Geometry\"   {.in2d } } \n",\
"	    {\"2D Geometry New\"   {.in2dnew } } \n",\
"	} \n",\
"\n",\
"	set ACISavailable [Ng_ACISCommand isACISavailable]\n",\
"	if {$ACISavailable == \"yes\" } {\n",\
"	    lappend types {\"ACIS Geometry\" {.sat} }\n",\
"	}\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	set file [tk_getOpenFile -filetypes $types]\n",\
"	if {$file != \"\"} {\n",\
"	    AddRecentFile $file\n",\
"	    Ng_LoadGeometry $file \n",\
"	    Ng_ParseGeometry\n",\
"	    set selectvisual geometry\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	    wm title . [concat \"$progname - \" $file]\n",\
"	    set dirname [file dirname $file]\n",\
"	    set basefilename [file tail [file rootname $file]]\n",\
"\n",\
"	    rebuildoccdialog\n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Save Geometry...\" \\\n",\
"    -command { \n",\
"	set occgeometryloaded [Ng_OCCCommand isoccgeometryloaded]\n",\
"	puts $occgeometryloaded\n",\
"	if {$occgeometryloaded == 1 } {\n",\
"	    set types {\n",\
"		{\"IGES Geometry file\"   {.igs} } \n",\
"		{\"STEP Geometry file\"   {.stp} } \n",\
"		{\"STL Geometry file\"   {.stl} } \n",\
"		{\"STL BIN Geometry file\"   {.stlb} } \n",\
"	    }\n",\
"	} {\n",\
"	    set types {\n",\
"		{\"STL Geometry file\"   {.stl} } \n",\
"		{\"STL BIN Geometry file\"   {.stlb} } \n",\
"	    }\n",\
"	}\n",\
"\n",\
"	set ACISavailable [Ng_ACISCommand isACISavailable]\n",\
"	puts $ACISavailable\n",\
"	if {$ACISavailable == \"yes\" } {\n",\
"	    lappend types {\"ACIS Geometry\" {.sat} }\n",\
"	}\n",\
"\n",\
"\n",\
"\n",\
"	set file [tk_getSaveFile -filetypes $types -initialdir $dirname -initialfile $basefilename ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_SaveGeometry $file \n",\
"	}\n",\
"    }\n",\
" \n",\
"\n",\
"\n",\
".ngmenu.file add cascade -label \"Recent Files\" -menu .ngmenu.file.recent \n",\
"menu .ngmenu.file.recent\n",\
"\n",\
"\n",\
"proc AddRecentFile { filename } {\n",\
"    global progname\n",\
"    global dirname\n",\
"    catch { [.ngmenu.file.recent delete $filename] }\n",\
"    .ngmenu.file.recent insert 0 command -label $filename \\\n",\
"	-command \"AddRecentFile {$filename}; \n",\
"                  Ng_LoadGeometry {$filename}; \n",\
"		  Ng_ParseGeometry;\n",\
"		  set selectvisual geometry;\n",\
"		  Ng_SetVisParameters;\n",\
"	          redraw;\n",\
"		  wm title . [concat \\\" $progname - $filename \\\"];\n",\
"                  set dirname {[file dirname $filename]};\n",\
"                  set basefilename {[file tail [file rootname $filename]]};\n",\
"                  rebuildoccdialog;\"\n",\
"    \n",\
"    if { [.ngmenu.file.recent index last] >= 6 } {\n",\
"	.ngmenu.file.recent delete last }\n",\
"    \n",\
"    saveinifile;\n",\
"    }\n",\
"loadinifile;\n",\
"\n",\
"\n",\
".ngmenu.file add separator\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Load Mesh...\" -accelerator \"<l><m>\" \\\n",\
"    -command {\n",\
"	set types {\n",\
"	    {\"Mesh file\"   {.vol}	} }\n",\
"	set file [tk_getOpenFile -filetypes $types -defaultextension \".vol\"]\n",\
"	if {$file != \"\"} {\n",\
"	    AddRecentMeshFile $file;\n",\
"	    Ng_LoadMesh $file; \n",\
"	    set selectvisual mesh\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	    Ng_ReadStatus; \n",\
"	    set dirname [file dirname $file]\n",\
"	    set basefilename [file tail [file rootname $file]]\n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add cascade -label \"Recent Meshes\" -menu .ngmenu.file.recentmesh \n",\
"menu .ngmenu.file.recentmesh\n",\
"\n",\
"\n",\
"proc AddRecentMeshFile { filename } {\n",\
"    global progname\n",\
"    global dirname\n",\
"    catch { [.ngmenu.file.recentmesh delete $filename] }\n",\
"    .ngmenu.file.recentmesh insert 0 command -label $filename \\\n",\
"	-command \"AddRecentMeshFile {$filename}; \n",\
"                  Ng_LoadMesh {$filename};\n",\
"		  set selectvisual mesh;\n",\
"		  Ng_SetVisParameters;\n",\
"	          redraw;\n",\
"		  wm title . [concat \\\" $progname - $filename \\\"];\n",\
"                  set dirname {[file dirname $filename]};\n",\
"                  set basefilename {[file tail [file rootname $filename]]};\n",\
"                  rebuildoccdialog;\"\n",\
"    \n",\
"    if { [.ngmenu.file.recentmesh index last] >= 6 } {\n",\
"	.ngmenu.file.recentmesh delete last }\n",\
"   \n",\
"    savemeshinifile;\n",\
"    }\n",\
"loadmeshinifile;\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Save Mesh...\" -accelerator \"<s><m>\" \\\n",\
"    -command {\n",\
"	set types {\n",\
"	    {\"Mesh file\"   {.vol}	} }\n",\
"\n",\
"	set file [tk_getSaveFile -filetypes $types -defaultextension \".vol\" -initialfile $basefilename -initialdir $dirname ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_SaveMesh $file }\n",\
"	AddRecentMeshFile $file;\n",\
"\n",\
"    }\n",\
"\n",\
".ngmenu.file add command -label \"Merge Mesh...\" \\\n",\
"    -command {\n",\
"	set types {\n",\
"	    {\"Mesh file\"   {.vol}	} }\n",\
"	set file [tk_getOpenFile -filetypes $types -defaultextension \".vol\"]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_MergeMesh $file; \n",\
"	    set selectvisual mesh\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	    Ng_ReadStatus; \n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Import Mesh...\" \\\n",\
"    -command { \n",\
"	set types {\n",\
"	    {\"Neutral format\"  {.mesh .emt} }\n",\
"	    {\"Surface mesh format\"  {.surf} }\n",\
"	    {\"Universal format\"  {.unv} }\n",\
"	    {\"Olaf format\"  {.emt} }\n",\
"	    {\"TET format\" {.tet} }\n",\
"	    {\"Pro/ENGINEER neutral format\" {.fnf} }\n",\
"	          }\n",\
"	set file [tk_getOpenFile -filetypes $types ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_ImportMesh $file \n",\
"	    set selectvisual mesh\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	    Ng_ReadStatus; \n",\
"	}\n",\
"    }\n",\
"\n",\
".ngmenu.file add command -label \"Export Mesh...\" \\\n",\
"    -command {\n",\
"	if { $exportfiletype == \"Elmer Format\" } {\n",\
"	    set file [tk_chooseDirectory]\n",\
"        } else {\n",\
"	    set file [tk_getSaveFile]\n",\
"	}\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_ExportMesh $file $exportfiletype \n",\
"	}\n",\
"    }\n",\
"\n",\
".ngmenu.file add cascade -label \"Export Filetype\" -menu .ngmenu.file.filetype \n",\
"\n",\
"menu .ngmenu.file.filetype \n",\
"\n",\
"\n",\
".ngmenu.file add separator\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Save Solution...\" \\\n",\
"    -command { \n",\
"	set types { \n",\
"            {\"Solution File\"  {.sol} } \n",\
"            {\"VTK File\"  {.vtk} } \n",\
"        }\n",\
"	set file [tk_getSaveFile -filetypes $types ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_SaveSolution $file \n",\
"	}\n",\
"    }\n",\
"\n",\
".ngmenu.file add command -label \"Import Solution...\" \\\n",\
"    -command { \n",\
"	set types { {\"Solution File\"  {.sol} } }\n",\
"	set file [tk_getOpenFile -filetypes $types -defaultextension \".sol\"  ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_ImportSolution $file \n",\
"	    set selectvisual solution\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set demostarttime [clock clicks -millisecond]\n",\
"set stopdemo 0\n",\
"proc demoredraw { } {\n",\
"    global demostarttime\n",\
"    global stopdemo\n",\
"    set curtime [clock clicks -millisecond]\n",\
"    set result [ Ng_DemoSetTime [expr $curtime - $demostarttime] ]\n",\
"    redraw\n",\
"    global videoactive\n",\
"    if { $videoactive == 1 } {\n",\
"        puts \"addframe\"\n",\
"        .ndraw Ng_VideoClip addframe\n",\
"    }\n",\
"    if { $result == 0 && $stopdemo == 0 } {\n",\
"	after 1 { demoredraw }\n",\
"    }\n",\
"}\n",\
".ngmenu.file add command -label \"Show Demo...\" \\\n",\
"    -command {\n",\
"	set types { {\"Demo File\"  {.dem} } }\n",\
"	set file [tk_getOpenFile -filetypes $types -defaultextension \".dem\"  ]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_ShowDemo $file \n",\
"	    set demostarttime [clock clicks -millisecond]\n",\
"	    set stopdemo 0\n",\
"	    demoredraw\n",\
" 	}\n",\
"     }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add separator\n",\
"\n",\
".ngmenu.file add command -label \"Snapshot...\" \\\n",\
"    -command { \n",\
"	set types { \n",\
"	    {\"JPG file\" {.jpg} } \n",\
"	    {\"GIF file\" {.gif} } \n",\
"	    {\"PPM file\" {.ppm} } \n",\
"	}\n",\
"	set file [tk_getSaveFile -filetypes $types]\n",\
"	if {$file != \"\"} {\n",\
"	    .ndraw Ng_SnapShot $file }\n",\
"    }\n",\
"\n",\
"\n",\
".ngmenu.file add cascade -label \"Video clip\" -menu .ngmenu.file.video\n",\
"menu .ngmenu.file.video\n",\
"\n",\
"set videoactive 0\n",\
".ngmenu.file.video add command -label \"start...\" \\\n",\
"    -command { \n",\
" 	set types { \n",\
" 	    {\"MPG file\" {.mpg} } \n",\
" 	}\n",\
" 	set file [tk_getSaveFile -filetypes $types]\n",\
" 	if {$file != \"\"} {\n",\
" 	    .ndraw Ng_VideoClip init $file \n",\
"            global videoactive\n",\
"            set videoactive 1\n",\
"        }\n",\
"     }\n",\
"\n",\
".ngmenu.file.video add command -label \"add frame...\" \\\n",\
"    -command {.ndraw Ng_VideoClip addframe }\n",\
"\n",\
".ngmenu.file.video add command -label \"one cycle\" \\\n",\
"    -command {\n",\
"	set visoptions.redrawperiodic 1\n",\
"	for { set j 0 } { $j < 100 } { incr j } {\n",\
"	    puts \"j =  $j\"\n",\
"	    Ng_Vis_Set time [expr (1000 * $j / 100)]\n",\
"	    redraw\n",\
"	    .ndraw Ng_VideoClip addframe \n",\
"	    after 200\n",\
"	}\n",\
"    }\n",\
"\n",\
".ngmenu.file.video add command -label \"finalize...\" \\\n",\
"    -command {\n",\
"        .ndraw Ng_VideoClip finalize \n",\
"        global videoactive\n",\
"        set videoactive 0\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Save Options\" \\\n",\
"    -command { saveoptions }\n",\
"\n",\
"\n",\
"    \n",\
"\n",\
".ngmenu.file add separator\n",\
"\n",\
"\n",\
".ngmenu.file add command -label \"Run tests ...\" \\\n",\
"    -command { runtestdialog }\n",\
".ngmenu.file add separator\n",\
"\n",\
".ngmenu.file add command -label \"Quit\" -accelerator \"<q>\" \\\n",\
"    -command { puts \"Thank you for using $progname\"; Ng_Exit; destroy . }\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.mesh\n",\
".ngmenu.mesh add command -label \"Generate Mesh\" -accelerator \"<g><m>\" \\\n",\
"    -command { \n",\
"	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}\n",\
"	Ng_ReadStatus\n",\
"	set selectvisual mesh\n",\
"	Ng_SetVisParameters\n",\
"	redraw\n",\
"    }\n",\
"\n",\
".ngmenu.mesh add command -label \"Stop Meshing\" \\\n",\
"    -command { Ng_StopMeshing }\n",\
"\n",\
".ngmenu.mesh add command -label \"Meshing Options...\" \\\n",\
"    -command meshingoptionsdialog\n",\
"\n",\
".ngmenu.mesh add separator\n",\
"\n",\
".ngmenu.mesh add command -label \"Delete Mesh\" \\\n",\
"    -command { Ng_New mesh; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.mesh add command -label \"Delete Vol Mesh\" \\\n",\
"    -command { Ng_DeleteVolMesh; Ng_ReadStatus; redraw }\n",\
"\n",\
"\n",\
".ngmenu.mesh add command -label \"Mesh Info\" \\\n",\
"    -command {\n",\
"	set dim [Ng_MeshInfo dim]\n",\
"	set np [Ng_MeshInfo np]\n",\
"	set ne [Ng_MeshInfo ne]\n",\
"	set nse [Ng_MeshInfo nse]\n",\
"	set nseg [Ng_MeshInfo nseg]\n",\
"	set bbox [Ng_MeshInfo bbox]\n",\
"	tk_messageBox -message  \"Dimension: $dim\\nPoints: $np\\nElements: $ne\\nSurface Els: $nse\\nSegments: $nseg\\nxmin [lindex $bbox 0] xmax [lindex $bbox 1]\\nymin [lindex $bbox 2] ymax [lindex $bbox 3]\\nzmin [lindex $bbox 4] zmax [lindex $bbox 5]\"\n",\
"    }\n",\
"\n",\
"\n",\
".ngmenu.mesh add command -label \"Mesh Quality\" \\\n",\
"    -command {\n",\
"	set inplanemin 0\n",\
"	set inplanemax 0\n",\
"	set betplanemin 0\n",\
"	set betplanemax 0\n",\
"	Ng_MeshQuality inplanemin inplanemax betplanemin betplanemax\n",\
"	puts \"Triangle angles : $inplanemin - $inplanemax\"\n",\
"	puts \"Tet angles      : $betplanemin - $betplanemax\"\n",\
"	tk_messageBox -message  \"Triangle angles : $inplanemin - $inplanemax \\n Tet angles      : $betplanemin - $betplanemax\"\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.mesh add command -label \"Check Surface Mesh\" \\\n",\
"    -command { Ng_CheckSurfaceMesh }\n",\
".ngmenu.mesh add command -label \"Check Volume Mesh\" \\\n",\
"    -command { Ng_CheckVolumeMesh }\n",\
"\n",\
".ngmenu.mesh add command -label \"Edit Boundary Conditions...\" \\\n",\
"    -command { bcpropdialog }\n",\
"\n",\
"if { $userlevel == 3 } {\n",\
"    .ngmenu.mesh add command -label \"Mesh Doctor...\" \\\n",\
"	-command { meshdoctordialog }\n",\
"}\n",\
"\n",\
".ngmenu.mesh add command -label \"METIS Mesh Partitioning...\" \\\n",\
"	-command { METISdialog }\n",\
"\n",\
".ngmenu.mesh add separator\n",\
"\n",\
".ngmenu.mesh add command -label \"Analyze Geometry\" \\\n",\
"    -command { Ng_GenerateMesh ag ag; Ng_ReadStatus; redraw }\n",\
".ngmenu.mesh add command -label \"Mesh Edges\" \\\n",\
"    -command { Ng_GenerateMesh me me; Ng_ReadStatus; redraw }\n",\
".ngmenu.mesh add command -label \"Mesh Surface\" \\\n",\
"    -command { set selectvisual mesh; Ng_SetVisParameters; \\\n",\
"		   Ng_GenerateMesh ms ms; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.mesh add command -label \"Optimize Surface\" \\\n",\
"    -command { Ng_GenerateMesh os os cmsmSm; redraw }\n",\
"\n",\
".ngmenu.mesh add cascade -label \"Surface Optim. Step\" -menu .ngmenu.mesh.surfoptstep \n",\
"\n",\
"menu .ngmenu.mesh.surfoptstep \n",\
".ngmenu.mesh.surfoptstep add command -label \"Mesh Smoothing\" \\\n",\
"    -command { Ng_GenerateMesh os os m; redraw}\n",\
".ngmenu.mesh.surfoptstep add command -label \"Edge swapping (topologic)\" \\\n",\
"    -command { Ng_GenerateMesh os os s; redraw}\n",\
".ngmenu.mesh.surfoptstep add command -label \"Edge swapping (metric)\" \\\n",\
"    -command { Ng_GenerateMesh os os S; redraw}\n",\
".ngmenu.mesh.surfoptstep add command -label \"Combine points\" \\\n",\
"    -command { Ng_GenerateMesh os os c; redraw}\n",\
"\n",\
"\n",\
".ngmenu.mesh add separator\n",\
".ngmenu.mesh add command -label \"Mesh Volume\" \\\n",\
"    -command { Ng_GenerateMesh mv mv; Ng_ReadStatus }\n",\
".ngmenu.mesh add command -label \"Optimize Volume\" \\\n",\
"    -command { Ng_GenerateMesh ov ov; Ng_ReadStatus }\n",\
".ngmenu.mesh add command -label \"Smooth Opt Volume\" \\\n",\
"    -command { Ng_GenerateMesh ov ov m; Ng_ReadStatus }\n",\
".ngmenu.mesh add command -label \"Smooth Opt Volume Jacobian\" \\\n",\
"    -command { Ng_GenerateMesh ov ov j; Ng_ReadStatus }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.geometry\n",\
".ngmenu.geometry add command -label \"Scan CSG Geometry\" -command { Ng_ParseGeometry }\n",\
".ngmenu.geometry add command -label \"CSG Options...\" -command geometryoptionsdialog\n",\
".ngmenu.geometry add command -label \"CSG Properties...\" \\\n",\
"    -command topleveldialog2 \n",\
"\n",\
".ngmenu.geometry add separator\n",\
"\n",\
".ngmenu.geometry add command -label \"STL Doctor...\" \\\n",\
"    -command { stldoctordialog; }\n",\
"\n",\
".ngmenu.geometry add command -label \"STL Info\" \\\n",\
"    -command {\n",\
"	set notriangles 0\n",\
"	set minx 0\n",\
"	set maxx 0\n",\
"	set miny 0\n",\
"	set maxy 0\n",\
"	set minz 0\n",\
"	set maxz 0\n",\
"	set trigscons 0\n",\
"	Ng_STLInfo notriangles minx maxx miny maxy minz maxz trigscons\n",\
"	set msgtext \"NO STL-Triangles : $notriangles\\nGeometry:\\nX = $minx - $maxx\\nY = $miny - $maxy\\nZ = $minz - $maxz\\nConsistency Check = $trigscons\\n\"\n",\
"	set msgtext \"$msgtext Status: [Ng_STLInfo status]\"\n",\
"	tk_messageBox -title \"STL Info\" -message  $msgtext -type ok \n",\
"    }\n",\
"\n",\
".ngmenu.geometry add separator\n",\
"\n",\
".ngmenu.geometry add command -label \"IGES/STEP Topology Explorer/Doctor...\" \\\n",\
"    -command { occdialog; }\n",\
"\n",\
"\n",\
".ngmenu.geometry add command -label \"OCC Construction\" \\\n",\
"    -command { Ng_OCCConstruction; }\n",\
"\n",\
"if { [Ng_ACISCommand isACISavailable] == \"yes\" } {\n",\
"    .ngmenu.geometry add command -label \"ACIS Topology Explorer...\" \\\n",\
"        -command { acisdialog; }\n",\
"\n",\
"    .ngmenu.geometry add command -label \"ACIS combine all\" \\\n",\
"        -command { Ng_ACISCommand combineall }\n",\
"    .ngmenu.geometry add command -label \"ACIS Create CT\" \\\n",\
"        -command { Ng_ACISCommand createct }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.view\n",\
".ngmenu.view add command -label \"Zoom all\" \\\n",\
"    -command { Ng_ZoomAll; redraw }\n",\
".ngmenu.view add command -label \"Center\" \\\n",\
"    -command { Ng_Center; redraw }\n",\
"\n",\
".ngmenu.view add command -label \"x-y plane\" \\\n",\
"    -command { Ng_StandardRotation xy; redraw }\n",\
".ngmenu.view add command -label \"y-x plane\" \\\n",\
"    -command { Ng_StandardRotation yx; redraw }\n",\
".ngmenu.view add command -label \"x-z plane\" \\\n",\
"    -command { Ng_StandardRotation xz; redraw }\n",\
".ngmenu.view add command -label \"z-x plane\" \\\n",\
"    -command { Ng_StandardRotation zx; redraw }\n",\
".ngmenu.view add command -label \"y-z plane\" \\\n",\
"    -command { Ng_StandardRotation yz; redraw }\n",\
".ngmenu.view add command -label \"z-y plane\" \\\n",\
"    -command { Ng_StandardRotation zy; redraw }\n",\
"\n",\
".ngmenu.view add command -label \"Viewing Options...\" \\\n",\
"    -command { viewingoptionsdialog; redraw }\n",\
".ngmenu.view add command -label \"Clipping Plane...\" \\\n",\
"    -command { clippingdialog; redraw }\n",\
".ngmenu.view add command -label \"Solution Data...\" \\\n",\
"    -command { visual_dialog; redraw }\n",\
".ngmenu.view add checkbutton -variable viewqualityplot \\\n",\
"    -label \"Quality Plot\" \\\n",\
"    -command { qualityviewdialog $viewqualityplot }\n",\
".ngmenu.view add checkbutton -variable memuseplot \\\n",\
"    -label \"Memory Usage\" \\\n",\
"    -command { memusedialog $memuseplot }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.meshsize\n",\
".ngmenu.meshsize add command -label \"Refine uniform\" \\\n",\
"    -command { Ng_Refine; Ng_HighOrder ${options.elementorder}; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.meshsize add command -label \"Second Order\" \\\n",\
"    -command { Ng_SecondOrder; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.meshsize add command -label \"Validate Second Order\" \\\n",\
"    -command { Ng_ValidateSecondOrder; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.meshsize add command -label \"High Order\" \\\n",\
"    -command { Ng_HighOrder ${options.elementorder}; Ng_ReadStatus; redraw }\n",\
"\n",\
".ngmenu.meshsize add separator\n",\
"\n",\
".ngmenu.meshsize add command -label \"Refinement Dialog...\" \\\n",\
"    -command { refinementdialog }\n",\
".ngmenu.meshsize add command -label \"Load Meshsize...\" \\\n",\
"    -command {\n",\
"	set types {\n",\
"	    {\"Meshsize file\"   {.msz}	} }\n",\
"	set file [tk_getOpenFile -filetypes $types]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_LoadMeshSize $file; \n",\
"	}\n",\
"    }\n",\
".ngmenu.meshsize add command -label \"MS from Surf Mesh\" \\\n",\
"    -command { Ng_MeshSizeFromSurfaceMesh }\n",\
"\n",\
"\n",\
"if { $userlevel == 3 } {\n",\
".ngmenu.meshsize add command -label \"Singular point ms\" \\\n",\
"    -command { Ng_SingularPointMS; }\n",\
"\n",\
".ngmenu.meshsize add command -label \"Singular edge ms\" \\\n",\
"    -command { Ng_SingularEdgeMS; }\n",\
"\n",\
".ngmenu.meshsize add separator\n",\
"\n",\
"set bisectfilename \"\";\n",\
"\n",\
".ngmenu.meshsize add command -label \"Bisection\" \\\n",\
"    -command { Ng_ReadStatus; set oldnp 0; set newnp $status_np; \n",\
"	Ng_ReadStatus;\n",\
"	\n",\
"	while { $oldnp < $newnp } {\n",\
"	    set level [expr $level+1]\n",\
"	    if { $bisectfilename == \"\"} {\n",\
"		Ng_Bisect;\n",\
"	    } else {\n",\
"		Ng_Bisect $bisectfilename;\n",\
"	    }\n",\
"	    Ng_ReadStatus;\n",\
"	    redraw; \n",\
"	    \n",\
"	    if { $bisectfilename == \"\"} {\n",\
"		set oldnp $newnp;\n",\
"		set newnp $status_np;\n",\
"		puts \"oldnp $oldnp newnp $newnp\";\n",\
"	    } else {\n",\
"		set oldnp $newnp;\n",\
"	    }\n",\
"	}\n",\
"    }\n",\
"\n",\
"}\n",\
"\n",\
".ngmenu.meshsize add command -label \"Load Refinement Info...\" \\\n",\
"    -command {\n",\
"	set types {\n",\
"	    {\"Refinement info\" {.refine} }}\n",\
"	set bisectfilename [tk_getOpenFile -filetypes $types]\n",\
"    }\n",\
"\n",\
".ngmenu.meshsize add command -label \"Z-Refinement\" \\\n",\
"    -command { Ng_ZRefinement 2; Ng_ReadStatus; redraw }\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.meshsize add cascade -label \"hp-Refinement\" -menu .ngmenu.meshsize.hpref\n",\
"menu .ngmenu.meshsize.hpref\n",\
".ngmenu.meshsize.hpref add command -label \"1 Level\" \\\n",\
"    -command { Ng_HPRefinement 1; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"2 Levels\" \\\n",\
"    -command { Ng_HPRefinement 2; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"3 Levels\" \\\n",\
"    -command { Ng_HPRefinement 3; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"4 Levels\" \\\n",\
"    -command { Ng_HPRefinement 4; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"5 Levels\" \\\n",\
"    -command { Ng_HPRefinement 5; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"6 Levels\" \\\n",\
"    -command { Ng_HPRefinement 6; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"7 Levels\" \\\n",\
"    -command { Ng_HPRefinement 7; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"8 Levels\" \\\n",\
"    -command { Ng_HPRefinement 8; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"9 Levels\" \\\n",\
"    -command { Ng_HPRefinement 9; Ng_ReadStatus; redraw }\n",\
".ngmenu.meshsize.hpref add command -label \"10 Levels\" \\\n",\
"    -command { Ng_HPRefinement 10; Ng_ReadStatus; redraw }\n",\
"\n",\
"\n",\
".ngmenu.meshsize add command -label \"Split to Tets\" \\\n",\
"    -command { Ng_Split2Tets; Ng_ReadStatus; redraw }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.special\n",\
".ngmenu.special add command -label \"Insert virtual boundary layer\" \\\n",\
"    -command { Ng_InsertVirtualBL; redraw }\n",\
".ngmenu.special add command -label \"Cut off and combine with other\" \\\n",\
"    -command { \n",\
"	set types { {\"Mesh file\"   {.vol}	} }\n",\
"	set file [tk_getOpenFile -filetypes $types]\n",\
"	if {$file != \"\"} {\n",\
"	    Ng_CutOffAndCombine $file;  }\n",\
"	redraw \n",\
"    }\n",\
".ngmenu.special add command -label \"Helmholtz Mesh grading\" \\\n",\
"    -command { Ng_HelmholtzMesh; }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"menu .ngmenu.help\n",\
".ngmenu.help add command -label \"Ng Help...\" \\\n",\
"	-command { help_main }\n",\
".ngmenu.view add checkbutton -label \"Help Line\" -variable showhelpline \\\n",\
"	-command {\n",\
"    if { $showhelpline == 1} {\n",\
"	pack .helpline -before .statbar -side bottom -fill x -padx 3p\n",\
"    } {\n",\
"	pack forget .helpline \n",\
"    }\n",\
"} \n",\
"\n",\
".ngmenu.help add command -label \"About...\" \\\n",\
"    -command {\n",\
"tk_messageBox -message \"This is NETGEN \\n mainly written by \\n Joachim Schöberl \\n thanks to \\n R. Gaisbauer, J. Gerstmayr\"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"frame .bubar -relief raised -bd 2\n",\
"pack .bubar -side top -fill x\n",\
"\n",\
"button .bubar.testb -text \"Test\" -command { Ng_SaveGeometry }\n",\
"button .bubar.surfm -text \"Generate Mesh\" -command \\\n",\
"    { set selectvisual mesh; \n",\
"	Ng_SetVisParameters;\n",\
"	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}\n",\
"	redraw \n",\
"    }\n",\
"button .bubar.stopm -text \"Stop\" -command \\\n",\
"    { Ng_StopMeshing;  set stopdemo 1 }\n",\
"button .bubar.exitb -text \"Quit\" \\\n",\
"    -command { .ngmenu.file invoke \"Quit\" }\n",\
"pack  .bubar.exitb .bubar.surfm .bubar.stopm -side left\n",\
"\n",\
"\n",\
"Ng_IsParallel;\n",\
"if { $parallel_netgen } {\n",\
"    button .bubar.visallb -text \"Parallel\" -command \\\n",\
"	{ paralleldialog; redraw } \n",\
"    pack .bubar.visallb -side left \n",\
"}\n",\
"\n",\
"button .bubar.zoomall -text \"Zoom All\" \\\n",\
"    -command { Ng_ZoomAll; redraw }\n",\
"\n",\
"button .bubar.center -text \"Center\" \\\n",\
"    -command { Ng_Center; redraw }\n",\
"\n",\
"tixOptionMenu .bubar.modesel \\\n",\
"    -options {\n",\
"	label.width  0\n",\
"	label.anchor e\n",\
"	menubutton.width 6\n",\
"    } \\\n",\
"    -variable drawmode\n",\
"\n",\
".bubar.modesel add command rotate -label Rotate\n",\
".bubar.modesel add command move -label Move\n",\
".bubar.modesel add command zoom -label Zoom\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set viewvals { geometry specpoints mesh solution}\n",\
"if { $userlevel == 3} {\n",\
"    set viewvals { geometry mesh specpoints surfmeshing modelview solution}\n",\
"}\n",\
"\n",\
"set viewvallabs(cross)     \"Cross\" \n",\
"set viewvallabs(geometry)  \"Geometry\" \n",\
"set viewvallabs(mesh)      \"Mesh\" \n",\
"set viewvallabs(specpoints) \"Edges\" \n",\
"set viewvallabs(surfmeshing) \"Mesh Gen\" \n",\
"set viewvallabs(modelview)     \"Modeller\" \n",\
"set viewvallabs(solution)     \"Solution\" \n",\
"\n",\
"tixOptionMenu .bubar.selview \\\n",\
"    -options {\n",\
"	label.width  0\n",\
"	label.anchor e\n",\
"	menubutton.width 10\n",\
"    } \\\n",\
"\n",\
"foreach viewv $viewvals {\n",\
"    .bubar.selview add command $viewv -label $viewvallabs($viewv)\n",\
"}\n",\
"\n",\
"\n",\
".bubar.selview config -variable selectvisual\n",\
".bubar.selview config -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"\n",\
"pack .bubar.modesel -side right\n",\
"pack forget .bubar.modesel\n",\
"pack .bubar.center .bubar.zoomall .bubar.selview -side right\n",\
"\n",\
".ngmenu.view add checkbutton -variable viewrotatebutton \\\n",\
"    -label \"Enable LeftButton Selection\" \\\n",\
"    -command { \n",\
"	if { $viewrotatebutton } {\n",\
"	    pack .bubar.modesel -side right\n",\
"	} {\n",\
"	    pack forget .bubar.modesel\n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"label .helpline -text \"None\"\n",\
"pack forget .helpline -side bottom -fill x\n",\
"\n",\
"frame .statbar -relief flat -bd 2\n",\
"pack .statbar -side bottom -fill x\n",\
"\n",\
"label .statbar.ptslabel -text \"Points: \"\n",\
"label .statbar.ptsval -textvariable status_np\n",\
"label .statbar.elslabel -text \"   Elements: \"\n",\
"label .statbar.elsval -textvariable status_ne\n",\
"label .statbar.selslabel -text \"   Surf Elements: \"\n",\
"label .statbar.selsval -textvariable status_nse\n",\
"label .statbar.task -textvariable status_task\n",\
"\n",\
"pack .statbar.ptslabel .statbar.ptsval -side left -ipady 3p \n",\
"pack .statbar.elslabel .statbar.elsval -side left -ipady 3p \n",\
"pack .statbar.selslabel .statbar.selsval -side left -ipady 3p\n",\
"\n",\
"\n",\
"\n",\
"tixMeter .statbar.per -value 0 -text 0%\n",\
".statbar.per configure -fillcolor blue\n",\
"\n",\
"pack .statbar.per -side right\n",\
"pack .statbar.task -side right -ipady 4\n",\
"\n",\
"set qualbaraxis(0) 0\n",\
"set qualbar(0) 0\n",\
"set qualbarnull(0) 0\n",\
"\n",\
"\n",\
"\n",\
"proc timer2 { } {\n",\
"    global status_np\n",\
"    global status_ne\n",\
"    global status_nse\n",\
"    global multithread_running\n",\
"    global multithread_redraw\n",\
"    global status_working\n",\
"    global status_task\n",\
"    global status_percent\n",\
"    global status_tetqualclasses\n",\
"    \n",\
"\n",\
"    Ng_ReadStatus \n",\
"\n",\
"    if { $multithread_redraw == 1 } {\n",\
"	set multithread_redraw 0;\n",\
"	redraw;\n",\
"        \n",\
"        global videoactive\n",\
"        if { $videoactive == 1 } {\n",\
"            puts \"addframe\"\n",\
"            .ndraw Ng_VideoClip addframe\n",\
"        }\n",\
"    }\n",\
"\n",\
"    global mem_moveable\n",\
"    set mem_moveable [Ng_MemInfo moveable]\n",\
"\n",\
"\n",\
"    .statbar.per config -value [expr $status_percent/100] -text [format %2.1f [expr 0.1*int(10*$status_percent)]]%\n",\
"\n",\
"\n",\
"    if { $multithread_running } {\n",\
"	pack .statbar.per -side right -before .statbar.task -padx 6\n",\
"    } { \n",\
"	pack forget .statbar.per\n",\
"    }\n",\
"	\n",\
"\n",\
"\n",\
"        if {[winfo exists .qualityview_dlg] == 1} {\n",\
"	\n",\
"	global qualbar\n",\
"	global qualbarnull\n",\
"	global qualbaraxis\n",\
"\n",\
"	set maxval 0\n",\
"	for {set i 0} {$i < 20} {incr i} {\n",\
"	    if {[lindex $status_tetqualclasses $i] > $maxval} {\n",\
"		set maxval [lindex $status_tetqualclasses $i]\n",\
"	    }\n",\
"	} \n",\
"\n",\
"	set ubound 1\n",\
"	while { $ubound < $maxval } {\n",\
"	    set ubound [expr {10 * $ubound}]\n",\
"	}\n",\
"	if { $ubound/5 > $maxval } {\n",\
"	    set ubound [expr $ubound/5]\n",\
"	}\n",\
"	if { $ubound/2 > $maxval } {\n",\
"	    set ubound [expr $ubound/2]\n",\
"	}\n",\
"\n",\
"\n",\
"	\n",\
"	for {set i 1} {$i <= 5} {incr i} {\n",\
"	    \n",\
"	    set value [expr { $i * $ubound / 5 }]\n",\
"	    .qualityview_dlg.c dchars $qualbaraxis($i) 0 end\n",\
"	    .qualityview_dlg.c insert $qualbaraxis($i) end $value  \n",\
"	}\n",\
"\n",\
"	\n",\
"	for {set i 0} {$i < 20} {incr i} {\n",\
"	    set x1 [expr {100 + ($i*15) + 2}]\n",\
"	    set x2 [expr {$x1+10}]\n",\
"	    \n",\
"	    set nbrs [lindex $status_tetqualclasses $i]\n",\
"	    set y [expr (249 - (200 * $nbrs / $ubound ) )]\n",\
"	    \n",\
"	    	    .qualityview_dlg.c coords $qualbar($i) $x1 250 $x2 $y\n",\
"\n",\
"	    if { $nbrs == 0 } {\n",\
"		.qualityview_dlg.c itemconfigure $qualbarnull($i) -text 0\n",\
"	    } {\n",\
"		.qualityview_dlg.c itemconfigure $qualbarnull($i) -text \"\" \n",\
"	    }		\n",\
"	}\n",\
"	\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"    if {[winfo exists .memuse_dlg] == 1} {    \n",\
"	\n",\
"	global memmark\n",\
"	set usemb [Ng_MemInfo usedmb]\n",\
"	for {set i 0} {$i < [string length $usemb] } {incr i} {\n",\
"	    if { [string index $usemb $i] == 0 } {\n",\
"		.memuse_dlg.c coords $memmark($i)  [expr 50+$i] 68 [expr 50+$i] 70\n",\
"	    } {\n",\
"		.memuse_dlg.c coords $memmark($i)  [expr 50+$i] 50 [expr 50+$i] 70\n",\
"	    }\n",\
"	}\n",\
"\n",\
"    }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    after 200 { timer2 }\n",\
"}\n",\
"timer2\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc bgerror { error } {\n",\
"    global errorInfo userlevel\n",\
"    if { $userlevel == 3} {\n",\
"	puts \"ERROR: $error\" \n",\
"	puts \"errinfo: $errorInfo\"\n",\
"    }\n",\
"    tk_messageBox -title \"Error Message\" -message $error -type ok \n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc smh2 { menuitem } {\n",\
"    if {[catch {$menuitem entrycget active -label} name]} {\n",\
"	set name \"    \"\n",\
"    } \n",\
"    show_menu_help $name \n",\
"    update idletasks\n",\
"}\n",\
"\n",\
"bind .ngmenu <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.file <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.geometry <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.mesh <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.view <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.meshsize <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.special <<MenuSelect>> { smh2 %W }\n",\
"bind .ngmenu.help <<MenuSelect>> { smh2 %W }\n",\
"\n",\
"\n",\
"bind . <q> { .ngmenu.file invoke \"Quit\" }\n",\
"bind . <l><g> { .ngmenu.file invoke \"Load Geometry...\" }  ; \n",\
"bind . <l><m> { .ngmenu.file invoke \"Load Mesh...\" }  ;\n",\
"bind . <s><m> { .ngmenu.file invoke \"Save Mesh...\" }  ;\n",\
"bind . <r><f> { .ngmenu.file activate \"Recent Files\" }  ;\n",\
"bind . <n><p> { newprimitivedialog }      ; bind . <e><p> { editprimitivedialog }\n",\
"bind . <e><s> { newsoliddialog }\n",\
"bind . <g><m> { .ngmenu.mesh invoke \"Generate Mesh\" }  ;\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc meshingoptionsdialog { } {\n",\
"\n",\
"    set w .options_dlg\n",\
"    \n",\
"    if {[winfo exists .options_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"\n",\
"\n",\
"	tixNoteBook $w.nb -ipadx 6 -ipady 6\n",\
"	\n",\
"	$w.nb add general -label \"General\" -underline 0\n",\
"	$w.nb add meshsize -label \"Mesh Size\"   -underline 0\n",\
"	$w.nb add chartopt -label \"STL Charts\" -underline 0\n",\
"	$w.nb add optimizer -label \"Optimizer\"   -underline 0\n",\
"	$w.nb add insider -label \"Insider\"   -underline 0\n",\
"	$w.nb add debug -label \"Debug\"   -underline 0\n",\
"\n",\
"\n",\
"	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	\n",\
"\n",\
"\n",\
"	\n",\
"	set f [$w.nb subwidget general]\n",\
"\n",\
"	set finevals { 1 2 3 4 5 6 }\n",\
"	set finelabs(1) \"very coarse\" \n",\
"	set finelabs(2) \"coarse\" \n",\
"	set finelabs(3) \"moderate\" \n",\
"	set finelabs(4) \"fine\" \n",\
"	set finelabs(5) \"very fine\" \n",\
"	set finelabs(6) \"user defined\"\n",\
"\n",\
"	tixOptionMenu $f.fine -label \"Mesh granularity : \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"\n",\
"\n",\
"	foreach finev $finevals {\n",\
"	    $f.fine add command $finev -label $finelabs($finev)\n",\
"	}\n",\
"	$f.fine config -variable meshoptions.fineness\n",\
"	$f.fine config -command { setgranularity }\n",\
"	global meshoptions.fineness\n",\
"	pack $f.fine\n",\
"\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	set mgsteps { ag me ms os mv ov }\n",\
"	set mgsteplabel(ag) \"Analyze Geometry\"\n",\
"	set mgsteplabel(me) \"Mesh Edges\"\n",\
"	set mgsteplabel(ms) \"Mesh Surface\"\n",\
"	set mgsteplabel(os) \"Optimize Surface\"\n",\
"	set mgsteplabel(mv) \"Mesh Volume\"\n",\
"	set mgsteplabel(ov) \"Optimize Volume\"\n",\
"\n",\
"	\n",\
"	tixOptionMenu $f.first -label \"First Step : \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"\n",\
"	tixOptionMenu $f.last -label \"Last Step : \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"\n",\
"	foreach step $mgsteps {\n",\
"	    $f.first add command $step -label $mgsteplabel($step)\n",\
"	    $f.last add command $step -label $mgsteplabel($step)\n",\
"	}\n",\
"\n",\
"	$f.first config -variable meshoptions.firststep \n",\
"	$f.last config  -variable meshoptions.laststep \n",\
"\n",\
"	pack $f.first $f.last\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	set msg(0) \"None\"\n",\
"	set msg(1) \"Least\"\n",\
"	set msg(2) \"Little\"\n",\
"	set msg(3) \"Moderate\"\n",\
"	set msg(4) \"Much\"\n",\
"	set msg(5) \"Most\"\n",\
"	\n",\
"	tixOptionMenu $f.msg -label \"Print Messages : \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"\n",\
"	foreach step {0 1 2 3 4 5 } {\n",\
"	    $f.msg add command $step -label $msg($step)\n",\
"	}\n",\
"\n",\
"	$f.msg config -variable options.printmsg \n",\
"	pack $f.msg\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	\n",\
"\n",\
"	checkbutton $f.parthread -text \"Parallel meshing thread\" \\\n",\
"	    -variable options.parthread\n",\
"	checkbutton $f.second -text \"Second order elements\" \\\n",\
"	    -variable options.secondorder\n",\
"	checkbutton $f.quad -text \"Quad dominated\" \\\n",\
"	    -variable options.quad -command {\n",\
"		if { ${options.quad} } {\n",\
"		    set meshoptions.laststep os\n",\
"		}\n",\
"	    }\n",\
"	checkbutton $f.invtets -text \"Invert volume elements\" \\\n",\
"	    -variable options.inverttets\n",\
"	checkbutton $f.invtrigs -text \"Invert surface elements\" \\\n",\
"	    -variable options.inverttrigs\n",\
"	checkbutton $f.azref -text \"Automatic Z-refinement\" \\\n",\
"	    -variable options.autozrefine\n",\
"\n",\
"	pack $f.parthread $f.second $f.quad $f.invtets $f.invtrigs $f.azref \n",\
"\n",\
"\n",\
"\n",\
"	tixControl $f.elementorder -label \"Element order: \" -integer true \\\n",\
"	    -variable options.elementorder -min 1 -max 20 \\\n",\
"	    -options {\n",\
"		entry.width 2\n",\
"		label.width 20\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"        pack $f.elementorder\n",\
"\n",\
"\n",\
"	tixControl $f.memory -label \"Large Memory \\[MB\\]: \" -integer true \\\n",\
"	    -variable options.memory -min 0 -max 2000 \\\n",\
"	    -options {\n",\
"		entry.width 5\n",\
"		label.width 20\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	global userlevel\n",\
"	if { $userlevel >= 3} { pack $f.memory }\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget meshsize]\n",\
"\n",\
"	tixControl $f.meshsize -label \"max mesh-size: \" -integer false \\\n",\
"	    -variable options.meshsize -min 1e-9 -max 1e6 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	tixControl $f.minmeshsize -label \"min mesh-size: \" -integer false \\\n",\
"	    -variable options.minmeshsize -min 0 -max 1e6 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	tixControl $f.grading -label \"mesh-size grading: \" -integer false \\\n",\
"	    -variable options.grading -min 0.1 -max 1 -step 0.1 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	\n",\
"	pack $f.meshsize $f.minmeshsize $f.grading\n",\
"\n",\
"\n",\
"\n",\
"	frame $f.msf \n",\
"	pack $f.msf\n",\
"\n",\
" 	tixLabelEntry $f.msf.ent -label \"mesh-size file: \"  \\\n",\
" 	    -labelside top \\\n",\
" 	    -options {  \n",\
" 		entry.textVariable options.meshsizefilename \n",\
" 		entry.width 35\n",\
" 		label.width 25\n",\
" 		label.anchor w\n",\
" 	    }	\n",\
" 	button $f.msf.btn -text \"Browse\" -command {\n",\
"	    global options.meshsizefilename\n",\
"	    set types {\n",\
"		{\"Meshsize file\"   {.msz}	} }\n",\
"	    set options.meshsizefilename [tk_getOpenFile -filetypes $types -initialfile ${options.meshsizefilename}]\n",\
"	}\n",\
"\n",\
" 	pack $f.msf.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4\n",\
" 	pack $f.msf.btn -side left -anchor s -padx 4 -pady 4\n",\
"\n",\
"\n",\
"	label $f.lab -text \"Additional mesh size restrictions:\"\n",\
"\n",\
"	\n",\
"	frame $f.csg -relief groove -borderwidth 3\n",\
"	pack $f.csg -fill x\n",\
"\n",\
"\n",\
"	frame $f.csg.curv\n",\
"	pack $f.csg.curv -anchor w\n",\
"\n",\
"	scale $f.csg.curv.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable options.curvaturesafety\n",\
"	label $f.csg.curv.la -text \"Elements per curvature radius\"\n",\
"	pack $f.csg.curv.sc $f.csg.curv.la -side left \n",\
"\n",\
"	frame $f.csg.elen\n",\
"	pack $f.csg.elen -anchor w\n",\
"	scale $f.csg.elen.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable options.segmentsperedge\n",\
"	label $f.csg.elen.la -text \"Elements per edge\"\n",\
"	pack $f.csg.elen.sc $f.csg.elen.la -side left\n",\
"	\n",\
"\n",\
"	\n",\
"	frame $f.stl -relief groove -borderwidth 3\n",\
"	pack $f.stl -fill x\n",\
"\n",\
"	frame $f.stl.r2\n",\
"	pack $f.stl.r2 -anchor w\n",\
"	scale $f.stl.r2.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthchartdistfac\n",\
"	checkbutton $f.stl.r2.bu -text \"STL - chart distance\" \\\n",\
"	    -variable stloptions.resthchartdistenable\n",\
"	pack $f.stl.r2.sc $f.stl.r2.bu -side left\n",\
"	\n",\
"	frame $f.stl.r6\n",\
"	pack $f.stl.r6 -anchor w\n",\
"	scale $f.stl.r6.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthlinelengthfac\n",\
"	checkbutton $f.stl.r6.bu -text \"STL - line length\" \\\n",\
"	    -variable stloptions.resthlinelengthenable\n",\
"	pack $f.stl.r6.sc $f.stl.r6.bu -side left\n",\
"	\n",\
"	frame $f.stl.r3\n",\
"	pack $f.stl.r3 -anchor w\n",\
"	scale $f.stl.r3.sc -orient horizontal -length 200 -from 0.2 -to 8 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthcloseedgefac\n",\
"	checkbutton $f.stl.r3.bu -text \"STL/IGES/STEP - close edges\" \\\n",\
"	-variable stloptions.resthcloseedgeenable\n",\
"	pack $f.stl.r3.sc $f.stl.r3.bu -side left\n",\
"	\n",\
"	frame $f.stl.r1\n",\
"	pack $f.stl.r1 -anchor w\n",\
"	scale $f.stl.r1.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthsurfcurvfac\n",\
"	checkbutton $f.stl.r1.bu -text \"STL - surface curvature\" \\\n",\
"	    -variable stloptions.resthsurfcurvenable\n",\
"	pack $f.stl.r1.sc $f.stl.r1.bu -side left\n",\
"\n",\
"	frame $f.stl.r3b\n",\
"	pack $f.stl.r3b -anchor w\n",\
"	scale $f.stl.r3b.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthedgeanglefac\n",\
"	checkbutton $f.stl.r3b.bu -text \"STL - edge angle\" \\\n",\
"	-variable stloptions.resthedgeangleenable\n",\
"	pack $f.stl.r3b.sc $f.stl.r3b.bu -side left\n",\
"	\n",\
"	frame $f.stl.r5\n",\
"	pack $f.stl.r5 -anchor w\n",\
"	scale $f.stl.r5.sc -orient horizontal -length 200 -from 0.2 -to 5 \\\n",\
"	    -resolution 0.1 -variable stloptions.resthsurfmeshcurvfac\n",\
"	checkbutton $f.stl.r5.bu -text \"STL - surface mesh curv\" \\\n",\
"	    -variable stloptions.resthsurfmeshcurvenable\n",\
"	pack $f.stl.r5.sc $f.stl.r5.bu -side left\n",\
"	\n",\
"	\n",\
"	checkbutton $f.stl.recalch -text \"STL - Recalc mesh size for surface optimization\" \\\n",\
"	    -variable stloptions.recalchopt\n",\
"	pack $f.stl.recalch\n",\
"\n",\
"	button $f.stl.calch -text \"Calc New H\" -command { redraw; Ng_STLCalcLocalH }\n",\
"	pack $f.stl.calch\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	set f [$w.nb subwidget chartopt]\n",\
"\n",\
"\n",\
"	label $f.lab1 -text \"Yellow Edges Angle ()\"\n",\
"	scale $f.scale1 -orient horizontal -length 300 \\\n",\
"	    -from 0 -to 90 -resolution 1  -tickinterval 10 \\\n",\
"	    -variable  stloptions.yangle \n",\
"\n",\
"	pack $f.lab1 $f.scale1\n",\
"\n",\
"	label $f.lab2e -text \"Edge Corner Angle ()\"\n",\
"	scale $f.scale2e -orient horizontal -length 360 -from 0 -to 180 \\\n",\
"	    -resolution 1  -tickinterval 20 \\\n",\
"	    -variable  stloptions.edgecornerangle \n",\
"	pack $f.lab2e $f.scale2e\n",\
"	\n",\
"	label $f.lab2 -text \"Chart Angle ()\"\n",\
"	scale $f.scale2 -orient horizontal -length 360 -from 0 -to 180 \\\n",\
"	    -resolution 1  -tickinterval 20 \\\n",\
"	    -variable  stloptions.chartangle \n",\
"	pack $f.lab2 $f.scale2\n",\
"	\n",\
"	label $f.lab2b -text \"Outer Chart Angle ()\"\n",\
"	scale $f.scale2b -orient horizontal -length 360 -from 0 -to 180 \\\n",\
"	    -resolution 1  -tickinterval 20 \\\n",\
"	    -variable  stloptions.outerchartangle \n",\
"	pack $f.lab2b $f.scale2b\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"		\n",\
"	set f [$w.nb subwidget optimizer]\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	tixControl $f.os2d -label \"Surface opt steps: \" -integer true \\\n",\
"	    -variable options.optsteps2d -min 0 -max 99 -step 1 \\\n",\
"	    -options {\n",\
"		entry.width 3\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	tixControl $f.os3d -label \"Volume opt steps: \" -integer true \\\n",\
"	    -variable options.optsteps3d -min 0 -max 99 -step 1 \\\n",\
"	    -options {\n",\
"		entry.width 3\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	\n",\
"	tixControl $f.elw -label \"Element size weight: \" -integer false \\\n",\
"	    -variable options.elsizeweight -min 0 -max 1 -step 0.1 \\\n",\
"	    -options {\n",\
"		entry.width 3\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	tixControl $f.wem -label \"Worst element measure: \" -integer false \\\n",\
"	    -variable options.opterrpow -min 1 -max 10 -step 1 \\\n",\
"	    -options {\n",\
"		entry.width 3\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	pack $f.os2d $f.os3d $f.elw $f.wem\n",\
"	\n",\
"	frame $f.badellimit\n",\
"	pack $f.badellimit -fill x\n",\
"	label $f.badellimit.lab -text \"bad element criterion\";\n",\
"	scale $f.badellimit.scale -orient horizontal -length 150 \\\n",\
"	    -from 160 -to 180 -resolution 1 \\\n",\
"	-variable options.badellimit\n",\
"	pack $f.badellimit.scale $f.badellimit.lab -side right -anchor s\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget insider]\n",\
"	\n",\
"\n",\
"\n",\
"	checkbutton $f.localh -text \"Use Local Meshsize\" \\\n",\
"	    -variable options.localh\n",\
"	checkbutton $f.delauney -text \"Use Delaunay\" \\\n",\
"	    -variable options.delaunay\n",\
"	checkbutton $f.checkoverlap -text \"Check Overlapping\" \\\n",\
"	    -variable options.checkoverlap\n",\
"	checkbutton $f.checkcb -text \"Check Chart Boundary\" \\\n",\
"	    -variable options.checkchartboundary\n",\
"	checkbutton $f.blockfill -text \"Do Blockfilling\" \\\n",\
"	    -variable options.blockfill\n",\
"\n",\
"	pack  $f.localh  $f.delauney $f.checkoverlap  $f.blockfill $f.checkcb   -anchor w\n",\
"\n",\
"\n",\
"\n",\
"	\n",\
"		set f [$w.nb subwidget debug]\n",\
"\n",\
"	frame $f.cb\n",\
"	pack $f.cb -side top\n",\
"\n",\
"	\n",\
"\n",\
"	checkbutton $f.cb.slowchecks -text \"Slow checks\" \\\n",\
"	    -variable debug.slowchecks -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.debugoutput -text \"Debugging outout\" \\\n",\
"	    -variable debug.debugoutput -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltexline -text \"Halt on exising line\" \\\n",\
"	    -variable debug.haltexistingline  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltoverlap -text \"Halt on Overlap\" \\\n",\
"	    -variable debug.haltoverlap  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltsuc -text \"Halt on success\" \\\n",\
"	    -variable debug.haltsuccess  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltnosuc -text \"Halt on no success\" \\\n",\
"	    -variable debug.haltnosuccess  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltlargequal -text \"Halt on large quality class\" \\\n",\
"	    -variable debug.haltlargequalclass  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltseg -text \"Halt on Segment:\" \\\n",\
"	    -variable debug.haltsegment  -command { Ng_SetDebugParameters }\n",\
"	checkbutton $f.cb.haltnode -text \"Halt on Node:\" \\\n",\
"	    -variable debug.haltnode  -command { Ng_SetDebugParameters }\n",\
"\n",\
"\n",\
"	pack $f.cb.slowchecks $f.cb.debugoutput $f.cb.haltexline $f.cb.haltoverlap $f.cb.haltsuc $f.cb.haltnosuc $f.cb.haltlargequal  $f.cb.haltseg   $f.cb.haltnode \n",\
"\n",\
"	frame $f.cb.hf\n",\
"	pack $f.cb.hf -pady 5\n",\
"	checkbutton $f.cb.hf.cb -text \"Halt on Face:\" \\\n",\
"	    -variable debug.haltface  -command { Ng_SetDebugParameters }\n",\
"	entry $f.cb.hf.ent -textvariable debug.haltfacenr -width 5  \n",\
"	pack $f.cb.hf.cb $f.cb.hf.ent -side left \n",\
"\n",\
"	checkbutton $f.cb.showactivechart -text \"Show Active Meshing-Chart\" \\\n",\
"	    -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	pack $f.cb.showactivechart\n",\
"	\n",\
"\n",\
"	frame $f.segs\n",\
"	pack $f.segs -pady 5\n",\
"	label $f.segs.lab1 -text \"P1:\";\n",\
"	entry $f.segs.ent1 -width 8 -relief sunken \\\n",\
"	    -textvariable debug.haltsegmentp1 \n",\
"	label $f.segs.lab2 -text \"P2:\";\n",\
"	entry $f.segs.ent2 -width 8 -relief sunken \\\n",\
"	    -textvariable debug.haltsegmentp2 \n",\
"	pack $f.segs.lab1 $f.segs.ent1 $f.segs.lab2 $f.segs.ent2  -side left\n",\
"\n",\
"\n",\
"\n",\
"	frame $f.cont -relief groove -borderwidth 3\n",\
"	pack $f.cont \n",\
"		\n",\
"	checkbutton $f.cont.multidrawing -text \"Draw Meshing\" \\\n",\
"	    -variable multithread_drawing \n",\
"	pack $f.cont.multidrawing\n",\
"	\n",\
"	checkbutton $f.cont.multitestmode -text \"Meshing Testmode\" \\\n",\
"	-variable multithread_testmode \n",\
"	pack $f.cont.multitestmode\n",\
"	\n",\
"	button $f.cont.goon -text \"Go On\" -command { set multithread_pause 0 }\n",\
"	pack $f.cont.multidrawing $f.cont.multitestmode $f.cont.goon -side left -expand yes\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"	global userlevel\n",\
"	if { $userlevel < 3} {\n",\
"	    $w.nb delete insider\n",\
"	    $w.nb delete debug\n",\
"	}\n",\
"\n",\
"\n",\
"\n",\
"		\n",\
"	\n",\
"		\n",\
"	\n",\
"\n",\
"	\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"	button $w.bu.apl -text \"Apply\" -command { \n",\
"	    [.options_dlg.nb subwidget meshsize].meshsize invoke\n",\
"	    [.options_dlg.nb subwidget meshsize].grading invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].os2d invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].os3d invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].elw invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].wem invoke\n",\
"\n",\
"	    Ng_SetMeshingParameters \n",\
"	    Ng_SetDebugParameters\n",\
"	}\n",\
"\n",\
"	button $w.bu.ok -text \"Done\" -command {\n",\
"	    [.options_dlg.nb subwidget meshsize].meshsize invoke\n",\
"	    [.options_dlg.nb subwidget meshsize].grading invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].os2d invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].os3d invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].elw invoke\n",\
"	    [.options_dlg.nb subwidget optimizer].wem invoke\n",\
"\n",\
"	    Ng_SetMeshingParameters\n",\
"	    Ng_SetDebugParameters\n",\
"	    wm withdraw .options_dlg\n",\
"	}\n",\
"\n",\
"	pack  $w.bu.apl $w.bu.ok -side left -expand yes\n",\
"    \n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Meshing Options\"\n",\
"	focus .options_dlg\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"meshingoptionsdialog\n",\
"wm withdraw .options_dlg\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc geometryoptionsdialog { } {\n",\
"\n",\
"\n",\
"    set w .geometry_dlg\n",\
"    \n",\
"    if {[winfo exists .geometry_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"	\n",\
"	global geooptions\n",\
"	\n",\
"	Ng_GeometryOptions get\n",\
"\n",\
"	checkbutton $w.drawcsg -text \"Draw Geometry\" \\\n",\
"	-variable geooptions.drawcsg \n",\
"	pack $w.drawcsg\n",\
"\n",\
"	frame $w.fac\n",\
"	pack $w.fac -pady 5\n",\
"	label $w.fac.lab -text \"Facets:\";\n",\
"	entry $w.fac.ent -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.facets\n",\
"	pack $w.fac.lab $w.fac.ent  -side left\n",\
"	\n",\
" \n",\
"	frame $w.det\n",\
"	pack $w.det -pady 5\n",\
"	label $w.det.lab -text \"Detail:\";\n",\
"	entry $w.det.ent -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.detail\n",\
"	pack $w.det.lab $w.det.ent  -side left\n",\
"	\n",\
"	frame $w.cox\n",\
"	pack $w.cox -pady 5\n",\
"	label $w.cox.lab -text \"min/max x:\";\n",\
"	entry $w.cox.ent1 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.minx\n",\
"	entry $w.cox.ent2 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.maxx\n",\
"	pack $w.cox.lab $w.cox.ent1 \\\n",\
"	    $w.cox.ent2  -side left\n",\
"	\n",\
"	frame $w.coy\n",\
"	pack $w.coy -pady 5\n",\
"	label $w.coy.lab -text \"min/max y:\";\n",\
"	entry $w.coy.ent1 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.miny\n",\
"	entry $w.coy.ent2 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.maxy\n",\
"	pack $w.coy.lab $w.coy.ent1 \\\n",\
"	    $w.coy.ent2  -side left\n",\
"	\n",\
"	frame $w.coz\n",\
"	pack $w.coz -pady 5\n",\
"	label $w.coz.lab -text \"min/max z:\";\n",\
"	entry $w.coz.ent1 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.minz\n",\
"	entry $w.coz.ent2 -width 8 -relief sunken \\\n",\
"	    -textvariable geooptions.maxz\n",\
"	pack $w.coz.lab $w.coz.ent1 \\\n",\
"	    $w.coz.ent2  -side left\n",\
"	\n",\
"\n",\
"\n",\
"	\n",\
"	\n",\
"\n",\
" 	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"\n",\
" 	button $w.bu.app -text \"Apply\" -command {\n",\
" 	    Ng_GeometryOptions set\n",\
" 	}\n",\
" 	button $w.bu.ok -text \"Done\" -command {\n",\
" 	    Ng_GeometryOptions set\n",\
" 	    destroy .geometry_dlg\n",\
" 	}\n",\
" 	pack  $w.bu.app $w.bu.ok -side left -expand yes\n",\
"    \n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Geometry options\"\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"proc viewingoptionsdialog { } {\n",\
"\n",\
"    global userlevel\n",\
"\n",\
"    set w .viewopts_dlg\n",\
"    \n",\
"    if {[winfo exists .viewopts_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
" \n",\
"\n",\
"\n",\
"\n",\
"	tixNoteBook $w.nb -ipadx 6 -ipady 6\n",\
"	\n",\
"	$w.nb add general -label \"General\" -underline 0\n",\
"	$w.nb add stl -label \"STL\" -underline 0\n",\
"	$w.nb add occ -label \"IGES/STEP\" -underline 0\n",\
"	$w.nb add mesh -label \"Mesh\"   -underline 0\n",\
"	$w.nb add light -label \"Light\"   -underline 0\n",\
"	$w.nb add edges -label \"Edges\"   -underline 0\n",\
"	$w.nb add misc -label \"Misc.\"  -underline 3\n",\
"\n",\
"\n",\
"	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget general]\n",\
"\n",\
"	checkbutton $f.backcol -text \"White Background\" \\\n",\
"	-variable viewoptions.whitebackground \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.cross -text \"Draw Coordinate Cross\" \\\n",\
"	-variable viewoptions.drawcoordinatecross \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.color -text \"Draw Color-bar\" \\\n",\
"	-variable viewoptions.drawcolorbar \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.netgen -text \"Draw Netgen-logo\" \\\n",\
"	-variable viewoptions.drawnetgenlogo \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	pack $f.backcol $f.cross $f.color $f.netgen\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget stl]\n",\
"\n",\
"	frame $f.show -relief groove -borderwidth 3\n",\
"	pack $f.show\n",\
"	checkbutton $f.show.showtrias -text \"Show STL-Triangles\" \\\n",\
"	    -variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }\n",\
"	pack $f.show.showtrias -anchor w\n",\
"	\n",\
"	checkbutton $f.show.showfilledtrias -text \"Show Filled Triangles\" \\\n",\
"	    -variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }\n",\
"	pack $f.show.showfilledtrias -anchor w\n",\
"	\n",\
"	checkbutton $f.show.showactivechart -text \"Show Active Meshing-Chart\" \\\n",\
"	    -variable stloptions.showactivechart -command { Ng_SetVisParameters; redraw }\n",\
"	pack $f.show.showactivechart -anchor w\n",\
"	\n",\
"	checkbutton $f.show.showedges -text \"Show Edges\" \\\n",\
"	    -variable stloptions.showedges -command { Ng_SetVisParameters; redraw }\n",\
"	pack $f.show.showedges -anchor w\n",\
"	\n",\
"	frame $f.special -relief groove -borderwidth 3\n",\
"	pack $f.special\n",\
"	checkbutton $f.special.showmarktrias -text \"Show Chart Triangles\" \\\n",\
"	    -variable stloptions.showmarktrias \\\n",\
"	    -command {set stldoctor.showfaces 0; Ng_STLDoctor; Ng_SetVisParameters; redraw }\n",\
"	pack $f.special.showmarktrias -side left\n",\
"\n",\
"	checkbutton $f.special.showfaces -text \"Show Faces\" \\\n",\
"	    -variable stldoctor.showfaces \\\n",\
"	    -command {set stloptions.showmarktrias 0; Ng_STLDoctor; Ng_SetVisParameters; redraw}    \n",\
"	pack $f.special.showfaces -side left\n",\
"\n",\
"	frame $f.fn -relief groove -borderwidth 3\n",\
"	pack $f.fn\n",\
"	label $f.fn.lab3 -text \"Chart/Face number:\"\n",\
"	scale $f.fn.scale3 -orient horizontal -length 200 -from 0 -to 200 \\\n",\
"	    -resolution 1  -tickinterval 50 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  stloptions.chartnumber \n",\
"	pack $f.fn.lab3 $f.fn.scale3 -side left\n",\
"	\n",\
"	frame $f.fo -relief groove -borderwidth 3\n",\
"	pack $f.fo\n",\
"	label $f.fo.lab -text \"Chart/Face Offset:\";\n",\
"	entry $f.fo.ent -width 5 -relief sunken \\\n",\
"	    -textvariable stloptions.chartnumberoffset\n",\
"	pack $f.fo.lab $f.fo.ent -side left\n",\
"\n",\
"	frame $f.mt\n",\
"	pack $f.mt -fill x\n",\
"	checkbutton $f.mt.bu -text \"Show Marked (Dirty) Triangles\" \\\n",\
"	    -variable stldoctor.showmarkedtrigs \\\n",\
"	    -command {Ng_STLDoctor; redraw}    \n",\
"	pack $f.mt.bu\n",\
"\n",\
"	frame $f.ep\n",\
"	pack $f.ep -fill x\n",\
"	checkbutton $f.ep.bu -text \"show edge corner points\" \\\n",\
"	    -variable stldoctor.showedgecornerpoints \\\n",\
"	    -command {Ng_STLDoctor; redraw}    \n",\
"	pack $f.ep.bu\n",\
"\n",\
"	frame $f.stt\n",\
"	pack $f.stt -fill x\n",\
"	checkbutton $f.stt.bu -text \"show touched triangle chart\" \\\n",\
"	    -variable stldoctor.showtouchedtrigchart \\\n",\
"	    -command {set stldoctor.showfaces 0; set stloptions.showmarktrias 1; \\\n",\
"			  Ng_STLDoctor; Ng_SetVisParameters; redraw}    \n",\
"	pack $f.stt.bu\n",\
"\n",\
"	frame $f.sml\n",\
"	pack $f.sml -fill x\n",\
"	checkbutton $f.sml.bu -text \"draw meshed edges\" \\\n",\
"	    -variable stldoctor.drawmeshededges \\\n",\
"	    -command {Ng_STLDoctor;}    \n",\
"	pack $f.sml.bu\n",\
"	\n",\
"	\n",\
"	frame $f.sm\n",\
"	pack $f.sm -fill x\n",\
"	checkbutton $f.sm.bu -text \"select with mouse\" \\\n",\
"	    -variable stldoctor.selectwithmouse\n",\
"	pack $f.sm.bu\n",\
"	\n",\
"	frame $f.st -relief groove -borderwidth 3\n",\
"	pack $f.st -fill x\n",\
"	label $f.st.lab -text \"Select triangle by number\";\n",\
"	entry $f.st.ent -width 5 -relief sunken \\\n",\
"	    -textvariable stldoctor.selecttrig\n",\
"	pack $f.st.ent $f.st.lab -side left -expand yes\n",\
"	\n",\
"	frame $f.vc -relief groove -borderwidth 3\n",\
"	pack $f.vc -fill x\n",\
"	checkbutton $f.vc.bu -text \"show vicinity\" \\\n",\
"	    -variable stldoctor.showvicinity \\\n",\
"	    -command {Ng_STLDoctor vicinity; redraw}\n",\
"	label $f.vc.lab -text \"vicinity size\";\n",\
"	scale $f.vc.sc -orient horizontal -length 200 -from 0 -to 200 \\\n",\
"	    -resolution 1 -variable stldoctor.vicinity \\\n",\
"	    -command { Ng_STLDoctor vicinity; redraw }\n",\
"	pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes\n",\
"	\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget occ]\n",\
"	\n",\
"	checkbutton $f.occshowsurfaces -text \"Show surfaces \" \\\n",\
"	    -variable occoptions.showsurfaces \\\n",\
"	    -command { Ng_SetOCCVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.occshowedges -text \"Show edges \" \\\n",\
"	    -variable occoptions.showedges \\\n",\
"	    -command { Ng_SetOCCVisParameters; redraw }\n",\
"\n",\
"	frame $f.deflection -relief groove -borderwidth 3\n",\
"	pack $f.deflection -fill x\n",\
"	button $f.deflection.lab -text \"Rebuild visualization data\" \\\n",\
"	    -command {\n",\
"		Ng_SetOCCVisParameters\n",\
"		Ng_OCCCommand buildvisualizationmesh\n",\
"		redraw\n",\
"	    }\n",\
"\n",\
"	tixControl $f.deflection.ent -label \"Visualization smoothness\" -integer false \\\n",\
"	    -variable occoptions.deflection -min 0.1 -max 3 -step 0.1 \\\n",\
"	    -options { entry.width 3 } \\\n",\
"	    -command { Ng_SetOCCVisParameters }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	pack $f.deflection.ent $f.deflection.lab -side left  -expand yes\n",\
"	pack $f.occshowsurfaces $f.occshowedges\n",\
"\n",\
"\n",\
"	\n",\
"	tixControl $f.showsolid -label \"Show solid (0 for all)\" -integer true \\\n",\
"            -variable occoptions.showsolidnr -min 0 -max 999 \\\n",\
"	    -options { entry.width 3 } \\\n",\
"	    -command { Ng_SetOCCVisParameters; redraw }\n",\
"    \n",\
"	tixControl $f.showsolid2 -label \"Show solid 2\" -integer true \\\n",\
"            -variable occoptions.showsolidnr2 -min 0 -max 999 \\\n",\
"	    -options { entry.width 3 } \\\n",\
"	    -command { Ng_SetOCCVisParameters; redraw }\n",\
"\n",\
"	button $f.subtract -text \"Subtract (2 minus 1)\" \\\n",\
"	    -command {\n",\
"		Ng_ACISCommand subtract ${occoptions.showsolidnr} ${occoptions.showsolidnr2}\n",\
"		redraw\n",\
"	    }\n",\
"\n",\
"	button $f.combine -text \"Combine all\" \\\n",\
"	    -command {\n",\
"		Ng_ACISCommand combineall\n",\
"		redraw\n",\
"	    }\n",\
"\n",\
"	pack $f.showsolid $f.showsolid2 $f.subtract $f.combine\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget mesh]\n",\
"\n",\
"	checkbutton $f.showcolor -text \"Colored Meshsize Visualization\" \\\n",\
"	    -variable viewoptions.colormeshsize \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"\n",\
"	checkbutton $f.showfilledtrigs -text \"Show filled triangles\" \\\n",\
"	-variable viewoptions.drawfilledtrigs \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"	checkbutton $f.showedges -text \"Show edges\" \\\n",\
"	-variable viewoptions.drawedges \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"\n",\
"	checkbutton $f.showoutline -text \"Show Triangle Outline\" \\\n",\
"	    -variable viewoptions.drawoutline \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	\n",\
"	tixControl $f.subdiv -label \"Subdivision\" -integer true \\\n",\
"            -variable visoptions.subdivisions -min 0 -max 8 \\\n",\
"	    -options { entry.width 2 } \\\n",\
"	    -command { puts \"mesh-subdivision\"; Ng_SetVisParameters; Ng_Vis_Set parameters; Ng_SetNextTimeStamp; redraw }\n",\
"	\n",\
"	\n",\
"	checkbutton $f.showbadels -text \"Show bad elements\" \\\n",\
"	    -variable viewoptions.drawbadels \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"	checkbutton $f.showprisms -text \"Show prisms\" \\\n",\
"	    -variable viewoptions.drawprisms \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.showpyramids -text \"Show pyramids\" \\\n",\
"	    -variable viewoptions.drawpyramids \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.showhexes -text \"Show hexes\" \\\n",\
"	    -variable viewoptions.drawhexes \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	frame $f.fshrink\n",\
"	label $f.fshrink.lab -text \"Shrink elements\"\n",\
"	scale $f.fshrink.scale -orient horizontal -length 200 -from 0 -to 1.0001 \\\n",\
"	    -resolution 0.01  -tickinterval 0.25 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.shrink\n",\
"\n",\
"\n",\
"	checkbutton $f.showidentified -text \"Show identified points\" \\\n",\
"	    -variable viewoptions.drawidentified \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.showmetispartition -text \"Show METIS Partition\" \\\n",\
"	    -variable viewoptions.drawmetispartition \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	checkbutton $f.showpointnumbers -text \"Show Point-numbers\" \\\n",\
"	    -variable viewoptions.drawpointnumbers \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showedgenumbers -text \"Show Edge-numbers\" \\\n",\
"	    -variable viewoptions.drawedgenumbers \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showfacenumbers -text \"Show Face-numbers\" \\\n",\
"	    -variable viewoptions.drawfacenumbers \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showelementnumbers -text \"Show Element-numbers\" \\\n",\
"	    -variable viewoptions.drawelementnumbers \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"	\n",\
"	tixControl $f.showdomain -label \"Show surface of domain\" -integer true \\\n",\
"            -variable viewoptions.drawdomainsurf -min 0 -max 50 \\\n",\
"	    -options { entry.width 2 } \\\n",\
"	    -command { Ng_SetVisParameters; Ng_Vis_Set parameters; redraw }\n",\
"    \n",\
"    \n",\
"    \n",\
"\n",\
"	frame $f.center -relief groove -borderwidth 3\n",\
"	pack $f.center -fill x\n",\
"	button $f.center.lab -text \"Set Center Point\" \\\n",\
"	    -command { Ng_SetVisParameters; Ng_Center; redraw }\n",\
"	entry $f.center.ent -width 5 -relief sunken \\\n",\
"	    -textvariable viewoptions.centerpoint \n",\
"	pack $f.center.ent $f.center.lab -side left  -expand yes\n",\
"	\n",\
"	frame $f.drawel -relief groove -borderwidth 3\n",\
"	pack $f.drawel -fill x\n",\
"	button $f.drawel.lab -text \"Draw Element\" \\\n",\
"	    -command { Ng_SetVisParameters; Ng_ZoomAll; redraw }\n",\
"	entry $f.drawel.ent -width 5 -relief sunken \\\n",\
"	    -textvariable viewoptions.drawelement \n",\
"	pack $f.drawel.ent $f.drawel.lab -side left  -expand yes\n",\
"\n",\
"        pack $f.showfilledtrigs\n",\
"	pack $f.showoutline $f.subdiv $f.showedges  $f.showbadels \n",\
"		pack $f.showdomain \n",\
"	pack $f.showpointnumbers \n",\
"	pack $f.showedgenumbers $f.showfacenumbers $f.showelementnumbers \n",\
"	pack $f.showmetispartition\n",\
"\n",\
"\n",\
"	frame $f.frametets\n",\
"	checkbutton $f.frametets.showtets -text \"Show Tets in domain \" \\\n",\
"	    -variable viewoptions.drawtets \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	tixControl $f.frametets.showtetsdomain -label \"\" -integer true \\\n",\
"	    -variable viewoptions.drawtetsdomain -min 0 -max 500 \\\n",\
"	    -options { entry.width 2 } \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	pack $f.frametets\n",\
"	pack $f.frametets.showtets $f.frametets.showtetsdomain -side left\n",\
"\n",\
"\n",\
"	pack $f.showcolor    $f.showpyramids $f.showprisms $f.showhexes $f.showidentified\n",\
"	\n",\
"	pack $f.fshrink \n",\
"	pack $f.fshrink.lab $f.fshrink.scale -side left\n",\
"	\n",\
"	\n",\
"	\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget light]\n",\
"	\n",\
"	label $f.lab1 -text \"Ambient Light\"\n",\
"	scale $f.scale1 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.amb \n",\
"	label $f.lab2 -text \"Diffuse Light\"\n",\
"	scale $f.scale2 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.diff \n",\
"	label $f.lab3 -text \"Specular Light\"\n",\
"	scale $f.scale3 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.light.spec \n",\
"	label $f.lab4 -text \"Material Shininess\"\n",\
"	scale $f.scale4 -orient horizontal -length 300 -from 0 -to 128 \\\n",\
"	    -resolution 1  -tickinterval 32 \\\n",\
"	    -command { Ng_SetVisParameters; redraw } -variable  viewoptions.mat.shininess \n",\
"	label $f.lab5 -text \"Material Transparency\"\n",\
"	scale $f.scale5 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	-command { Ng_SetVisParameters; redraw } -variable  viewoptions.mat.transp \n",\
"	\n",\
"	pack $f.lab1 $f.scale1 $f.lab2 $f.scale2 $f.lab3 $f.scale3 $f.lab4 $f.scale4 $f.lab5 $f.scale5\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget edges]\n",\
"\n",\
"	checkbutton $f.showedges -text \"Show Edges\" \\\n",\
"	    -variable viewoptions.drawededges \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showpoints -text \"Show Points\" \\\n",\
"	    -variable viewoptions.drawedpoints \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showpointnrs -text \"Show Points Nrs\" \\\n",\
"	    -variable viewoptions.drawedpointnrs \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.showtang -text \"Show CP Tangents\" \\\n",\
"	    -variable viewoptions.drawedtangents \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	checkbutton $f.drawedgenrs -text \"Show Edge Nrs\" \\\n",\
"	    -variable viewoptions.drawededgenrs \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"	pack $f.showedges $f.showpoints $f.showpointnrs $f.showtang $f.drawedgenrs\n",\
"\n",\
"	frame $f.center -relief groove -borderwidth 3\n",\
"	pack $f.center -fill x\n",\
"	button $f.center.lab -text \"Set Center Point\" \\\n",\
"	    -command { Ng_SetVisParameters; Ng_Center; redraw }\n",\
"	entry $f.center.ent -width 5 -relief sunken \\\n",\
"	    -textvariable viewoptions.centerpoint \n",\
"	pack $f.center.ent $f.center.lab -side left  -expand yes\n",\
"	\n",\
"\n",\
"\n",\
"	frame $f.f1\n",\
"	pack $f.f1 -pady 5\n",\
"	label $f.f1.lab -text \"SpecPoint Veclen\"\n",\
"	entry $f.f1.ent -width 5 -relief sunken -textvariable viewoptions.specpointvlen\n",\
"	pack $f.f1.lab $f.f1.ent\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"		set f [$w.nb subwidget misc]\n",\
"\n",\
"	frame $f.point -relief groove -borderwidth 3\n",\
"\n",\
"	frame $f.point.dp\n",\
"	\n",\
"	checkbutton $f.point.dp.drawpoint -text \"Draw Point\" \\\n",\
"	    -variable viewoptions.drawspecpoint \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	entry $f.point.dp.px -width 8 -relief sunken -textvariable viewoptions.specpointx\n",\
"	entry $f.point.dp.py -width 8 -relief sunken -textvariable viewoptions.specpointy\n",\
"	entry $f.point.dp.pz -width 8 -relief sunken -textvariable viewoptions.specpointz\n",\
"\n",\
"	pack $f.point.dp.drawpoint $f.point.dp.px $f.point.dp.py $f.point.dp.pz -side left\n",\
"\n",\
"	pack $f.point.dp\n",\
"\n",\
"	checkbutton $f.point.center -text \"Use as Center\" \\\n",\
"	    -variable viewoptions.usecentercoords \\\n",\
"	    -command { \n",\
"		if { ${viewoptions.usecentercoords} } {\n",\
"		    set viewoptions.centerx ${viewoptions.specpointx}\n",\
"		    set viewoptions.centery ${viewoptions.specpointy}\n",\
"		    set viewoptions.centerz ${viewoptions.specpointz}\n",\
"		    Ng_SetVisParameters; Ng_Center\n",\
"		    redraw\n",\
"		} {\n",\
"		    Ng_SetVisParameters\n",\
"		}\n",\
"		\n",\
"		    \n",\
"	    }\n",\
"\n",\
"	pack $f.point.center\n",\
"	\n",\
"	pack $f.point -fill x -ipady 3\n",\
"\n",\
"\n",\
"	\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"\n",\
"	button $w.bu.done -text \"Done\" -command {\n",\
"	    Ng_SetVisParameters;\n",\
"	    redraw\n",\
"	    destroy .viewopts_dlg\n",\
"	}\n",\
"	button $w.bu.apply -text \"Apply\" -command {\n",\
"	    Ng_SetVisParameters;\n",\
"	    redraw\n",\
"	}\n",\
"	pack $w.bu.apply $w.bu.done -expand yes -side left\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Viewing options\"\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set clippingdialog_pop1 0\n",\
"set clippingdialog_pop2 0\n",\
"set clippingdialog_pop3 0\n",\
"set clippingdialog_pop4 0\n",\
"\n",\
"\n",\
"proc clippingdialog { } {\n",\
"\n",\
"    global clippingdialog_pop1\n",\
"    global clippingdialog_pop2\n",\
"    global clippingdialog_pop3\n",\
"    global clippingdialog_pop4\n",\
"    set clippingdialog_pop1 1\n",\
"    set clippingdialog_pop2 1\n",\
"    set clippingdialog_pop3 1\n",\
"    set clippingdialog_pop4 1\n",\
"    \n",\
"    set w .clipping_dlg\n",\
"    \n",\
"    if {[winfo exists .clipping_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"\n",\
"	label $w.lab1 -text \"Normal x\"\n",\
"	scale $w.scale1 -orient horizontal -length 300 -from -1 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.5 \\\n",\
"	    -variable  viewoptions.clipping.nx \\\n",\
"	    -command { popupcheckredraw2 clippingdialog_pop1 ${viewoptions.clipping.enable} }\n",\
"\n",\
"	\n",\
"	label $w.lab2 -text \"Normal y\"\n",\
"	scale $w.scale2 -orient horizontal -length 300 -from -1 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.5 \\\n",\
"	    -variable  viewoptions.clipping.ny \\\n",\
"	    -command { popupcheckredraw2 clippingdialog_pop2 ${viewoptions.clipping.enable} }\n",\
"\n",\
"	label $w.lab3 -text \"Normal z\"\n",\
"	scale $w.scale3 -orient horizontal -length 300 -from -1 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.5 \\\n",\
"	    -variable  viewoptions.clipping.nz \\\n",\
"	    -command { popupcheckredraw2 clippingdialog_pop3 ${viewoptions.clipping.enable} }\n",\
"	label $w.lab4 -text \"Distance\"\n",\
"	scale $w.scale4 -orient horizontal -length 300 -from -1 -to 1.001 \\\n",\
"	    -resolution 0.0001  -tickinterval 0.5 \\\n",\
"	    -variable  viewoptions.clipping.dist \\\n",\
"	    -command { popupcheckredraw2 clippingdialog_pop4 ${viewoptions.clipping.enable} }\n",\
"	\n",\
"	\n",\
"	tixControl $w.clipdomain -label \"Clip only domain\" -integer true \\\n",\
"	    -variable viewoptions.clipping.onlydomain -min 0 -max 50 \\\n",\
"	    -options { entry.width 2 } \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"	tixControl $w.donotclipdomain -label \"Do not clip domain\" -integer true \\\n",\
"	    -variable viewoptions.clipping.notdomain -min 0 -max 50 \\\n",\
"	    -options { entry.width 2 } \\\n",\
"	    -command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	pack $w.lab1 $w.scale1 $w.lab2 $w.scale2 $w.lab3 $w.scale3 $w.lab4 $w.scale4 $w.clipdomain $w.donotclipdomain\n",\
"\n",\
"	\n",\
"	checkbutton $w.cb1 -text \"Enable clipping\" \\\n",\
"	    -variable viewoptions.clipping.enable \\\n",\
"	    -command { Ng_SetVisParameters; redraw } \n",\
"	\n",\
"	pack $w.cb1\n",\
"	\n",\
"\n",\
"	\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"	button $w.bu.cancle -text \"Done\" -command \"destroy $w\"\n",\
"	pack $w.bu.cancle  -expand yes\n",\
"	\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Clipping Plane\"\n",\
"		focus $w\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc refinementdialog { } {\n",\
"\n",\
"    set w .refinement_dlg\n",\
"    \n",\
"    if {[winfo exists .refinement_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"\n",\
"	\n",\
"	tixControl $w.meshsize -label \"max mesh-size: \" -integer false \\\n",\
"	    -variable options.meshsize -min 1e-6 -max 1e6 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	pack $w.meshsize\n",\
"\n",\
"	global localh\n",\
"	set localh 1\n",\
"	tixControl $w.loch -label \"local mesh-size: \" -integer false \\\n",\
"	    -variable localh -min 1e-6 -max 1e6 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	\n",\
"	pack $w.loch\n",\
"	\n",\
"	\n",\
"	button $w.restface -text \"Restrict H at face\"  \\\n",\
"	    -command {\n",\
"		.refinement_dlg.meshsize invoke\n",\
"		.refinement_dlg.loch invoke\n",\
"		Ng_RestrictH face $localh\n",\
"	    }\n",\
"	button $w.restedge -text \"Restrict H at edge\"  \\\n",\
"	    -command {\n",\
"		.refinement_dlg.meshsize invoke\n",\
"		.refinement_dlg.loch invoke\n",\
"		Ng_RestrictH edge $localh\n",\
"	    }\n",\
"	button $w.restelement -text \"Restrict H at element\"  \\\n",\
"	    -command {\n",\
"		.refinement_dlg.meshsize invoke\n",\
"		.refinement_dlg.loch invoke\n",\
"		Ng_RestrictH element $localh\n",\
"	    }\n",\
"	button $w.restpoint -text \"Restrict H at point\"  \\\n",\
"	    -command {\n",\
"		.refinement_dlg.meshsize invoke\n",\
"		.refinement_dlg.loch invoke\n",\
"		Ng_RestrictH point $localh\n",\
"	    }\n",\
"\n",\
"\n",\
"	pack $w.restface $w.restedge $w.restelement $w.restpoint\n",\
"\n",\
"\n",\
"\n",\
"	button $w.anisoedge -text \"Declare Anisotropic edge\"  \\\n",\
"	    -command {\n",\
"		Ng_Anisotropy edge \n",\
"	    }\n",\
"	pack $w.anisoedge\n",\
"	\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"\n",\
"	button $w.bu.cancle -text \"Done\" -command \"destroy .refinement_dlg\"\n",\
"	button $w.bu.refine -text \"Refine\"  \\\n",\
"	    -command { \n",\
"		set oldnp 0; set newnp $status_np; \n",\
"		while { $oldnp < $newnp } {\n",\
"		    set level [expr $level+1]\n",\
"		    Ng_Bisect; \n",\
"		    Ng_HighOrder ${options.elementorder}\n",\
"		    Ng_ReadStatus;\n",\
"		redraw; \n",\
"		    set oldnp $newnp\n",\
"		    set newnp $status_np\n",\
"		    puts \"oldnp $oldnp newnp $newnp\"\n",\
"		}\n",\
"	    }	 \n",\
"	button $w.bu.zrefine -text \"Z-Refine\"  \\\n",\
"	    -command { Ng_ZRefinement; Ng_ReadStatus; redraw; }\n",\
"   \n",\
"	pack $w.bu.zrefine $w.bu.refine $w.bu.cancle  -expand yes -side left\n",\
"		\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Select Refinement\"\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc bcpropdialog { } {\n",\
"\n",\
"    set w .bcprop_dlg\n",\
"    \n",\
"    if {[winfo exists .bcprop_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	frame $w.face  -borderwidth 3\n",\
"	pack $w.face -fill x\n",\
"	label $w.face.lab -text \"face index:\"\n",\
"	label $w.face.ent -text 1 -padx 4\n",\
"	button $w.face.next -text \"next\" -command {\n",\
"	    set w .bcprop_dlg;	\n",\
"	    set facenr [$w.face.ent cget -text]\n",\
"	    if {$facenr == [Ng_BCProp getnfd]} {\n",\
"		set facenr 1 \n",\
"	    } {\n",\
"		set facenr [expr $facenr + 1]\n",\
"	    }\n",\
"	    $w.face.ent configure -text $facenr\n",\
"	    Ng_BCProp setactive $facenr\n",\
"	    set bcnr [Ng_BCProp getbc $facenr]\n",\
"	    $w.bc.ent delete 0 end\n",\
"	    $w.bc.ent insert 0 $bcnr\n",\
"\n",\
"	    redraw\n",\
"	} \n",\
"	button $w.face.prev -text \"prev\" -command {\n",\
"	    set w .bcprop_dlg;	\n",\
"	    set facenr [$w.face.ent cget -text]\n",\
"	    if {$facenr == 1} {\n",\
"		set facenr [Ng_BCProp getnfd]\n",\
"	    } {\n",\
"		set facenr [expr $facenr - 1]\n",\
"	    }\n",\
"	    $w.face.ent configure -text $facenr\n",\
"	    Ng_BCProp setactive $facenr\n",\
"	    set bcnr [Ng_BCProp getbc $facenr]\n",\
"	    $w.bc.ent delete 0 end\n",\
"	    $w.bc.ent insert 0 $bcnr\n",\
"\n",\
"	    redraw\n",\
"	} \n",\
"	\n",\
"	\n",\
"	pack $w.face.lab $w.face.ent $w.face.prev $w.face.next  -side left  \n",\
"	\n",\
"	frame $w.bc  -borderwidth 3\n",\
"	pack $w.bc -fill x\n",\
"	label $w.bc.lab -text \"bc property:\"\n",\
"	entry $w.bc.ent -width 5 -relief sunken \n",\
"	button $w.bc.but -text \"change\" -command { \n",\
"	    set w .bcprop_dlg;	    \n",\
"	    Ng_BCProp setbc [$w.face.ent cget -text] [$w.bc.ent get]; \n",\
"	}\n",\
"	button $w.bc.but2 -text \"all\" -command { \n",\
"	    set w .bcprop_dlg;	    \n",\
"	    Ng_BCProp setall [$w.bc.ent get]; \n",\
"	}\n",\
"	pack $w.bc.lab $w.bc.ent $w.bc.but $w.bc.but2 -side left  -expand yes\n",\
"\n",\
"	frame $w.bcname  -borderwidth 3\n",\
"	pack $w.bcname -fill x\n",\
"	label $w.bcname.lab -text \"bc name:\"\n",\
"	label $w.bcname.ent -text \"-\"\n",\
"	pack $w.bcname.lab $w.bcname.ent -side left  -expand yes\n",\
"	\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"	button $w.bu.close -text \"Close\" -command { destroy .bcprop_dlg }\n",\
"\n",\
"	pack $w.bu.close  -expand yes -side left\n",\
"		\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Boundary Conditions\"\n",\
"    }\n",\
"\n",\
"    focus $w \n",\
"\n",\
"    set facenr [Ng_BCProp getactive]\n",\
"    $w.face.ent configure -text $facenr\n",\
"    \n",\
"    set bcnr [Ng_BCProp getbc $facenr]\n",\
"    $w.bc.ent delete 0 end\n",\
"    $w.bc.ent insert 0 $bcnr\n",\
"\n",\
"    set bcname [Ng_BCProp getbcname $facenr]\n",\
"    $w.bcname.ent configure -text $bcname\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc METISdialog { } {\n",\
"\n",\
"    set w .metis_dlg\n",\
"    set w.parts 64\n",\
"    \n",\
"    if {[winfo exists .metis_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	frame $w.a -borderwidth 0\n",\
"	frame $w.b -borderwidth 0\n",\
"	pack $w.a $w.b\n",\
"\n",\
"	label $w.a.lab -text \"Number of partitions:\"\n",\
"	entry $w.a.ent -textvariable w.parts -width 4 -relief sunken\n",\
"\n",\
"	button $w.b.start -text \"Start METIS\" -command { \n",\
"	    Ng_Metis ${w.parts}\n",\
"	    redraw\n",\
"	}\n",\
"	button $w.b.cancel -text \"Cancel\" -command { destroy .metis_dlg }\n",\
"	pack $w.a.lab $w.a.ent -side left  -expand yes\n",\
"	pack $w.b.start $w.b.cancel -side left\n",\
"\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"METIS Partitioning\"\n",\
"	focus $w\n",\
" \n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc stloptionsdialog { } {\n",\
"\n",\
"    set w .stlopts_dlg\n",\
"    \n",\
"    if {[winfo exists .stlopts_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	tixNoteBook $w.nb -ipadx 6 -ipady 6\n",\
"			\n",\
"						\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	\n",\
"	\n",\
"	\n",\
"	\n",\
"\n",\
"	\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x -ipady 3\n",\
"\n",\
"\n",\
"    button $w.bu.apply -text \"Apply\" -command { redraw; Ng_GenerateMesh 1 2}\n",\
"	button $w.bu.cancle -text \"Done\" -command { destroy .stlopts_dlg }\n",\
"    pack $w.bu.cancle  $w.bu.apply  -side left -expand yes\n",\
"    \n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"	wm title $w \"STL Options\"\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"proc stldoctordialog { } {\n",\
"\n",\
"    set wd .stldoctor_dlg\n",\
"\n",\
"    if {[winfo exists .stldoctor_dlg] == 1} {\n",\
"	wm withdraw $wd\n",\
"	wm deiconify $wd\n",\
"	focus $wd \n",\
"    } {\n",\
"	\n",\
"    toplevel $wd\n",\
"\n",\
"    tixNoteBook $wd.nb -ipadx 6 -ipady 6\n",\
"	\n",\
"    $wd.nb add general -label \"General\" -underline 0\n",\
"    $wd.nb add topology -label \"Edit Topology\"  -underline 5\n",\
"    $wd.nb add edges -label \"Edit Edges\"   -underline 5\n",\
"    $wd.nb add normals -label \"Edit Normals\"   -underline 5\n",\
"    $wd.nb add advanced -label \"Advanced\"   -underline 0\n",\
"\n",\
"\n",\
"    pack $wd.nb -expand yes -fill both -padx 5 -pady 5 -side top	\n",\
"\n",\
"\n",\
"    \n",\
"    set f [$wd.nb subwidget general]\n",\
"\n",\
"\n",\
"    frame $f.show\n",\
"    pack $f.show -fill x\n",\
"    checkbutton $f.show.showtrias -text \"Show STL-Triangles\" \\\n",\
"	-variable stloptions.showtrias -command { Ng_SetVisParameters; redraw }\n",\
"    pack $f.show.showtrias -anchor w\n",\
"    \n",\
"    checkbutton $f.show.showfilledtrias -text \"Show Filled Triangles\" \\\n",\
"	-variable stloptions.showfilledtrias -command { Ng_SetVisParameters; redraw }\n",\
"    pack $f.show.showfilledtrias -anchor w\n",\
"\n",\
"    set selmodevals { 0 1 2 3 4 }\n",\
"    set selmodelabs(0) \"triangle\" \n",\
"    set selmodelabs(1) \"edge\" \n",\
"    set selmodelabs(2) \"point\" \n",\
"    set selmodelabs(3) \"line\" \n",\
"    set selmodelabs(4) \"line cluster\" \n",\
"\n",\
"    tixOptionMenu $f.selmode -label \"Double Click selects :\" \\\n",\
"	-options {\n",\
"	    label.width  19\n",\
"	    label.anchor e\n",\
"	    menubutton.width 15\n",\
"	} \n",\
"\n",\
"    foreach selmodev $selmodevals {\n",\
"	$f.selmode add command $selmodev -label $selmodelabs($selmodev)\n",\
"    }\n",\
"    $f.selmode config -variable stldoctor.selectmode\n",\
"    $f.selmode config -command { Ng_STLDoctor }\n",\
"    global stldoctor.selectmode\n",\
"    pack $f.selmode\n",\
"\n",\
"    frame $f.sm\n",\
"    pack $f.sm -fill x\n",\
"    checkbutton $f.sm.bu -text \"select with mouse\" \\\n",\
"	-variable stldoctor.selectwithmouse\n",\
"    pack $f.sm.bu \n",\
"\n",\
"    frame $f.st -relief groove -borderwidth 3\n",\
"    pack $f.st -fill x\n",\
"    label $f.st.lab -text \"Select triangle by number\";\n",\
"    entry $f.st.ent -width 5 -relief sunken \\\n",\
"	-textvariable stldoctor.selecttrig\n",\
"    pack $f.st.ent $f.st.lab -side left -expand yes\n",\
"\n",\
"    frame $f.vc -relief groove -borderwidth 3\n",\
"    pack $f.vc -fill x\n",\
"    checkbutton $f.vc.bu -text \"show vicinity\" \\\n",\
"	-variable stldoctor.showvicinity \\\n",\
"	-command {Ng_STLDoctor vicinity; redraw}\n",\
"    label $f.vc.lab -text \"vicinity size\";\n",\
"    scale $f.vc.sc -orient horizontal -length 200 -from 0 -to 200 \\\n",\
"	-resolution 1 -variable stldoctor.vicinity \\\n",\
"	-command { Ng_STLDoctor vicinity; redraw }\n",\
"    pack $f.vc.bu $f.vc.lab $f.vc.sc -expand yes\n",\
"\n",\
"    frame $f.ge -relief groove -borderwidth 3\n",\
"    pack $f.ge -fill x\n",\
"    button $f.ge.neighbourangles -text \"calc neighbourangles\" -command {Ng_STLDoctor neighbourangles}\n",\
"    button $f.ge.showcoords -text \"show coords of touched triangle\" -command {Ng_STLDoctor showcoords}\n",\
"    button $f.ge.moveptm -text \"move point to middle of trianglepoints\" -command {Ng_STLDoctor movepointtomiddle; redraw}\n",\
"    button $f.ge.destroy0trigs -text \"destroy 0-volume triangles\" -command {Ng_STLDoctor destroy0trigs}\n",\
"    pack $f.ge.neighbourangles $f.ge.showcoords $f.ge.moveptm $f.ge.destroy0trigs -expand yes \n",\
"\n",\
"\n",\
"    button $f.ge.cancle -text \"Done\" -command {destroy .stldoctor_dlg }\n",\
"    pack $f.ge.cancle -expand yes\n",\
"\n",\
"        set f [$wd.nb subwidget topology]\n",\
"\n",\
"    frame $f.oc -relief groove -borderwidth 3\n",\
"    pack $f.oc -fill x\n",\
"    button $f.oc.bu -text \"invert orientation of selected trig\" -command {Ng_STLDoctor invertselectedtrig; redraw }\n",\
"    button $f.oc.bu2 -text \"orient after selected trig\" -command {Ng_STLDoctor orientafterselectedtrig; redraw }\n",\
"    pack $f.oc.bu $f.oc.bu2 -side left  -expand yes\n",\
"\n",\
"    button $f.toperr -text \"mark inconsistent triangles\" -command {Ng_STLDoctor marktoperrortrigs; redraw }\n",\
"\n",\
"    button $f.deltrig -text \"delete selected triangle\" -command {Ng_STLDoctor deleteselectedtrig; redraw }\n",\
"    button $f.geosmooth -text \"geometric smoothing\" -command {Ng_STLDoctor smoothgeometry; redraw }\n",\
"\n",\
"    pack $f.toperr $f.deltrig $f.geosmooth\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"        set f [$wd.nb subwidget edges]\n",\
"\n",\
"\n",\
"    frame $f.be -relief groove -borderwidth 3 \n",\
"    pack $f.be -fill x\n",\
"    label $f.be.lab -text \"build edges with yellow angle:\";\n",\
"    scale $f.be.sc -orient horizontal -length 200 -from 0 -to 100 \\\n",\
"	-resolution 0.5\n",\
"    $f.be.sc config -variable stloptions.yangle \n",\
"    $f.be.sc config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }\n",\
"    label $f.be.lab2 -text \"continue edges with yellow angle:\";\n",\
"    scale $f.be.sc2 -orient horizontal -length 200 -from 0 -to 100 \\\n",\
"	-resolution 0.5\n",\
"    $f.be.sc2 config -variable stloptions.contyangle \n",\
"    $f.be.sc2 config -command { Ng_SetSTLParameters; Ng_STLDoctor buildedges; redraw }\n",\
"\n",\
"\n",\
"\n",\
"    button $f.be.buildedges -text \"Build Edges\" -command {Ng_STLDoctor buildedges; redraw}\n",\
"    pack $f.be.lab $f.be.sc $f.be.lab2 $f.be.sc2 $f.be.buildedges -expand yes\n",\
"\n",\
"    frame $f.se\n",\
"    pack $f.se -fill x\n",\
"    checkbutton $f.se.bu -text \"show excluded\" \\\n",\
"	-variable stldoctor.showexcluded \\\n",\
"	-command {Ng_STLDoctor; redraw}\n",\
"    pack $f.se.bu \n",\
"\n",\
"    \n",\
"    set edgeselmodevals { 0 1 2 3 4 }\n",\
"    set edgeselmodelabs(0) \"no change\" \n",\
"    set edgeselmodelabs(1) \"undefined\" \n",\
"    set edgeselmodelabs(2) \"confirmed\" \n",\
"    set edgeselmodelabs(3) \"candidate\"\n",\
"    set edgeselmodelabs(4) \"excluded\"\n",\
"\n",\
"    tixOptionMenu $f.edgeselmode -label \"Double Click sets edge :\" \\\n",\
"	-options {\n",\
"	    label.width  19\n",\
"	    label.anchor e\n",\
"	    menubutton.width 15\n",\
"	} \n",\
"\n",\
"    foreach edgeselmodev $edgeselmodevals {\n",\
"	$f.edgeselmode add command $edgeselmodev -label $edgeselmodelabs($edgeselmodev)\n",\
"    }\n",\
"    $f.edgeselmode config -variable stldoctor.edgeselectmode\n",\
"    $f.edgeselmode config -command { Ng_STLDoctor }\n",\
"    global stldoctor.edgeselectmode\n",\
"    pack $f.edgeselmode\n",\
"\n",\
"    \n",\
"    frame $f.edg -relief groove -borderwidth 3\n",\
"    pack $f.edg -fill x\n",\
"\n",\
"\n",\
"\n",\
"    frame $f.edg.f0\n",\
"    pack $f.edg.f0\n",\
"    button $f.edg.f0.confirmedge -text \"confirm\" -command {Ng_STLDoctor confirmedge; redraw}\n",\
"    button $f.edg.f0.candidateedge -text \"candidate\" -command {Ng_STLDoctor candidateedge; redraw}\n",\
"    button $f.edg.f0.excludeedge -text \"exclude\" -command {Ng_STLDoctor excludeedge; redraw}\n",\
"    button $f.edg.f0.undefinededge -text \"undefined\" -command {Ng_STLDoctor undefinededge; redraw}\n",\
"    pack $f.edg.f0.confirmedge $f.edg.f0.candidateedge $f.edg.f0.excludeedge $f.edg.f0.undefinededge  -side left\n",\
"\n",\
"    frame $f.edg.fa\n",\
"    pack $f.edg.fa\n",\
"    button $f.edg.fa.setallundefined -text \"all undefined\" -command {Ng_STLDoctor setallundefinededges; redraw}\n",\
"    button $f.edg.fa.erasecandidates -text \"candidates to undefined\" -command {Ng_STLDoctor erasecandidateedges; redraw}\n",\
"    pack $f.edg.fa.setallundefined $f.edg.fa.erasecandidates -side left\n",\
"\n",\
"\n",\
"    frame $f.edg.fb\n",\
"    pack $f.edg.fb\n",\
"    button $f.edg.fb.confirmcandidates -text \"candidates to confirmed\" -command {Ng_STLDoctor confirmcandidateedges; redraw}\n",\
"    button $f.edg.fb.confirmedtocandidates -text \"confirmed to candidates\" -command {Ng_STLDoctor confirmedtocandidateedges; redraw}\n",\
"    pack $f.edg.fb.confirmcandidates $f.edg.fb.confirmedtocandidates -side left\n",\
"\n",\
"    frame $f.edg.f1\n",\
"    frame $f.edg.f2\n",\
"    frame $f.edg.f3\n",\
"    frame $f.edg.f4\n",\
"    pack $f.edg.f1 $f.edg.f2 $f.edg.f3 $f.edg.f4\n",\
"\n",\
"    button $f.edg.f1.exportedges -text \"export edges\" -command {Ng_STLDoctor exportedges}\n",\
"    button $f.edg.f1.importedges -text \"import edges\" -command {Ng_STLDoctor importedges; redraw}\n",\
"    button $f.edg.f1.saveedgedata -text \"save edgedata\" \\\n",\
"	-command { \n",\
"	    set types {\n",\
"		{\"Netgen Edgedata\"   {.ned} } \n",\
"	    }\n",\
"	    set file [tk_getSaveFile -filetypes $types -defaultextension \".ned\"]\n",\
"	    if {$file != \"\"} {\n",\
"		Ng_STLDoctor saveedgedata $file\n",\
"	}\n",\
"    }\n",\
"\n",\
"    button $f.edg.f1.loadedgedata -text \"load edgedata\" \\\n",\
"	-command { \n",\
"	    set types {\n",\
"		{\"Netgen Edgedata\"  {.ned} }\n",\
"	    }\n",\
"	    set file [tk_getOpenFile -filetypes $types -defaultextension \".ned\"]\n",\
"	    if {$file != \"\"} {\n",\
"		Ng_STLDoctor loadedgedata $file \n",\
"		puts \"loading done\"\n",\
"		\n",\
"		redraw\n",\
"		\n",\
"	    }\n",\
"	} \n",\
"\n",\
"    button $f.edg.f1.importAVLedges -text \"import AVL edges\" \\\n",\
"	-command {\n",\
"	    set types {{\"Edge file\"  {.edg }}}\n",\
"\n",\
"	    set file [tk_getOpenFile -filetypes $types -defaultextension \".edg\"]\n",\
"	    if {$file != \"\"} {\n",\
"		Ng_STLDoctor importexternaledges $file; \n",\
"	    }\n",\
"	}\n",\
"\n",\
"    pack $f.edg.f1.importAVLedges $f.edg.f1.loadedgedata $f.edg.f1.saveedgedata -side left\n",\
"\n",\
"    frame $f.edg2 -relief groove -borderwidth 3\n",\
"    pack $f.edg2 -fill x\n",\
"\n",\
"\n",\
"    label $f.edg2.lab -text \"length (%):\"\n",\
"    scale $f.edg2.sc -orient horizontal -length 200 -from 0 -to 100 \\\n",\
"	-resolution 0.5 \\\n",\
"        -variable stldoctor.longlinefact \n",\
"\n",\
"     button $f.edg2.undoedge -text \"undo last edge change\" -command {Ng_STLDoctor undoedgechange; redraw}\n",\
"    \n",\
"      pack $f.edg2.undoedge -expand yes\n",\
"\n",\
"\n",\
"\n",\
"        set f [$wd.nb subwidget normals]\n",\
"\n",\
"    frame $f.dt -relief groove -borderwidth 3\n",\
"    pack $f.dt -fill x\n",\
"    label $f.dt.lab -text \"dirty triangle factor\";\n",\
"    entry $f.dt.ent -width 5 -relief sunken \\\n",\
"	-textvariable stldoctor.dirtytrigfact\n",\
"    pack $f.dt.ent $f.dt.lab -side left  -expand yes\n",\
"\n",\
"    frame $f.srt -relief groove -borderwidth 3\n",\
"    pack $f.srt -fill x\n",\
"    button $f.srt.bu -text \"smooth reverted triangles geometric\" -command {Ng_STLDoctor smoothrevertedtrigs; redraw }\n",\
"    entry $f.srt.ent -width 5 -relief sunken \\\n",\
"	-textvariable stldoctor.smoothangle\n",\
"    pack $f.srt.ent $f.srt.bu -side left  -expand yes\n",\
"\n",\
"    frame $f.bdt -relief groove -borderwidth 3\n",\
"    pack $f.bdt -fill x\n",\
"    button $f.bdt.bu -text \"mark dirty triangles\" -command {Ng_STLDoctor markdirtytrigs; redraw }\n",\
"    button $f.bdt.bu2 -text \"smooth dirty triangles normal\" -command {Ng_STLDoctor smoothdirtytrigs; redraw }\n",\
"    pack $f.bdt.bu $f.bdt.bu2 -side left  -expand yes\n",\
"\n",\
"    \n",\
"    frame $f.sno -relief groove -borderwidth 3\n",\
"    pack $f.sno\n",\
"    \n",\
"    label $f.sno.labrough -text \"rough\"\n",\
"    scale $f.sno.scsmooth -orient horizontal -length 100 -from 0 -to 0.8 \\\n",\
"	-resolution 0.01 -variable stldoctor.smoothnormalsweight \\\n",\
"	-command { Ng_SetSTLParameters }\n",\
"    label $f.sno.labsmooth -text \"smooth\"\n",\
"    button $f.sno.smoothnormals -text \"smooth normals\" -command { Ng_STLDoctor smoothnormals; redraw}\n",\
"\n",\
"\n",\
"\n",\
"    pack $f.sno.labrough $f.sno.scsmooth $f.sno.labsmooth $f.sno.smoothnormals -side left -padx 5\n",\
"\n",\
"    frame $f.no -relief groove -borderwidth 3\n",\
"    pack $f.no -fill x\n",\
"\n",\
"    button $f.no.marknonsmoothnormals -text \"mark non-smooth triangles\" -command {Ng_STLDoctor marknonsmoothnormals; redraw}\n",\
"    button $f.no.calcnormals -text \"calculate normals from geometry\" -command {Ng_STLDoctor calcnormals; redraw}\n",\
"\n",\
"    pack  $f.no.marknonsmoothnormals $f.no.calcnormals -expand yes\n",\
"\n",\
"\n",\
"        set f [$wd.nb subwidget advanced]\n",\
"\n",\
"\n",\
"    frame $f.sc\n",\
"    pack $f.sc -fill x\n",\
"    checkbutton $f.sc.bu -text \"spiral check\" \\\n",\
"	-variable stldoctor.spiralcheck \\\n",\
"	-command {Ng_STLDoctor;}    \n",\
"    checkbutton $f.sc.bu2 -text \"cone check\" \\\n",\
"	-variable stldoctor.conecheck \\\n",\
"	-command {Ng_STLDoctor;}    \n",\
"    pack $f.sc.bu $f.sc.bu2\n",\
"\n",\
"\n",\
"    tixControl $f.gtol -label \"load-geometry tolerance factor\" -integer false \\\n",\
"	-variable stldoctor.geom_tol_fact \\\n",\
"	-options {\n",\
"	    entry.width 8\n",\
"	    label.width 30\n",\
"	    label.anchor e\n",\
"	}	\n",\
"    pack $f.gtol\n",\
"\n",\
"    button $f.adap -text \"Apply\" -command {\n",\
"	[.stldoctor_dlg.nb subwidget advanced].gtol invoke\n",\
"	Ng_STLDoctor; \n",\
"    }\n",\
"    pack $f.adap -expand yes\n",\
"\n",\
"\n",\
"        wm withdraw $wd\n",\
"    wm geom $wd +100+100\n",\
"    wm deiconify $wd\n",\
"    wm title $wd \"STL Doctor\"\n",\
"\n",\
"    focus $wd\n",\
"}\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc meshdoctordialog { } {\n",\
"\n",\
"    set w .meshdoc_dlg\n",\
"    global meshdoctor.active\n",\
"\n",\
"    if {[winfo exists .meshdoc_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	set meshdoctor.active 1\n",\
"	Ng_MeshDoctor;\n",\
"\n",\
"\n",\
"	frame $w.vis -relief groove -borderwidth 3\n",\
"	pack $w.vis\n",\
"\n",\
"	checkbutton $w.vis.showfilledtrigs -text \"Show filled triangles\" \\\n",\
"	-variable viewoptions.drawfilledtrigs \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"	checkbutton $w.vis.showedges -text \"Show edges\" \\\n",\
"	-variable viewoptions.drawedges \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"	\n",\
"\n",\
"	checkbutton $w.vis.showoutline -text \"Show Triangle Outline\" \\\n",\
"	-variable viewoptions.drawoutline \\\n",\
"	-command { Ng_SetVisParameters; redraw }\n",\
"\n",\
"	pack $w.vis.showfilledtrigs  $w.vis.showoutline $w.vis.showedges\n",\
"\n",\
"	tixControl $w.markedgedist -label \"Mark edge dist: \" -integer true \\\n",\
"	    -min 0 -max 999  \\\n",\
"	    -variable meshdoc.markedgedist \\\n",\
"	    -options {\n",\
"		entry.width 3\n",\
"		label.width 20\n",\
"		label.anchor e\n",\
"	    } \\\n",\
"	    -command {\n",\
"		Ng_MeshDoctor markedgedist ${meshdoc.markedgedist}\n",\
"		redraw\n",\
"	    }\n",\
"	pack $w.markedgedist\n",\
"	\n",\
"	button $w.deledge -text \"Delete marked segments\" -command {\n",\
"	    Ng_MeshDoctor deletemarkedsegments\n",\
"	    redraw\n",\
"	}\n",\
"	pack $w.deledge\n",\
"	\n",\
"	button $w.close -text \"Close\" -command { \n",\
"	    set meshdoctor.active 0;\n",\
"	    Ng_MeshDoctor;\n",\
"	    destroy .meshdoc_dlg \n",\
"	}\n",\
"	pack $w.close -expand yes\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Mesh Doctor\"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc qualityviewdialog { show } {\n",\
"\n",\
"    set w .qualityview_dlg\n",\
"    \n",\
"    if {[winfo exists .qualityview_dlg] == 1} {\n",\
"\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw .qualityview_dlg\n",\
"	    wm deiconify $w\n",\
"	    focus $w \n",\
"	} {\n",\
"	    wm withdraw $w\n",\
"	}\n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	set c $w.c\n",\
"\n",\
"	canvas $c -relief raised -width 450 -height 300\n",\
"	pack $w.c -side top -fill x\n",\
"\n",\
"	set plotFont {Helvetica 12}\n",\
"	set smallFont {Helvetica 12}\n",\
"\n",\
"	$c create line 100 250 400 250 -width 2\n",\
"	$c create line 100 250 100 50 -width 2\n",\
"\n",\
"	for {set i 0} {$i <= 10} {incr i} {\n",\
"	    set x [expr {100 + ($i*30)}]\n",\
"	    $c create line $x 250 $x 245 -width 2\n",\
"	    if { [expr {$i % 2}] == 0 } {\n",\
"		$c create text $x 254 -text [format %1.1f [expr 0.1*$i]] -anchor n -font $plotFont\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	global qualbar\n",\
"	global qualbarnull\n",\
"	global qualbaraxis\n",\
"\n",\
"	for {set i 0} {$i <= 5} {incr i} {\n",\
"	    set y [expr {250 - ($i*40)}]\n",\
"	    $c create line 100 $y 105 $y -width 2\n",\
"\n",\
"	    set qualbaraxis($i) \\\n",\
"		[$c create text 96 $y -text [expr $i*50].0 -anchor e -font $plotFont]\n",\
"	}\n",\
"\n",\
"	for {set i 0} {$i < 20} {incr i} {\n",\
"	    set x1 [expr {100 + ($i*15) + 2}]\n",\
"	    set x2 [expr {$x1+10}]\n",\
"	    set y [expr {250 - 10 * $i}]\n",\
"	    set qualbar($i) [$c create rectangle $x1 250 $x2 245 -fill blue]\n",\
"	    set qualbarnull($i) [$c create text [expr {($x1+$x2)/2}] 245 -text 0 -anchor s -font $smallFont -fill blue]	\n",\
"	}\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu\n",\
"		\n",\
"	button $w.close -text \"Close\" \\\n",\
"	    -command { \n",\
"		wm withdraw .qualityview_dlg\n",\
"		set viewqualityplot 0\n",\
"	    }\n",\
"	pack $w.close\n",\
"    \n",\
"	\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw $w\n",\
"	    wm geom $w +100+100\n",\
"	    wm deiconify $w\n",\
"	    wm title $w \"Mesh Quality\"\n",\
"	    focus $w\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc memusedialog { show } {\n",\
"\n",\
"    set w .memuse_dlg\n",\
"    \n",\
"    if {[winfo exists .memuse_dlg] == 1} {\n",\
"\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw .memuse_dlg\n",\
"	    wm deiconify $w\n",\
"	    focus $w \n",\
"	} {\n",\
"	    wm withdraw $w\n",\
"	}\n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	set c $w.c\n",\
"\n",\
"	canvas $c -relief raised -width 600 -height 300\n",\
"	pack $w.c -side top -fill x\n",\
"\n",\
"	set plotFont {Helvetica 18}\n",\
"	set smallFont {Helvetica 12}\n",\
"\n",\
"\n",\
"	global memmark\n",\
"	for {set i 0} {$i < 512} { incr i } {\n",\
"	    set memmark($i) [$c create line [expr 50+$i] 50 [expr 50+$i] 70 -fill blue]\n",\
"	}\n",\
"\n",\
"\n",\
"	set plotFont {Helvetica 18}\n",\
"	set smallFont {Helvetica 12}\n",\
"\n",\
"	$c create text 50 90 -text \"0 GB\" -anchor n -font $plotFont\n",\
"	$c create text 178 90 -text \"1 GB\" -anchor n -font $plotFont\n",\
"	$c create text 306 90 -text \"2 GB\" -anchor n -font $plotFont\n",\
"	$c create text 434 90 -text \"3 GB\" -anchor n -font $plotFont\n",\
"	$c create text 562 90 -text \"4 GB\" -anchor n -font $plotFont\n",\
"\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu\n",\
"		\n",\
"	button $w.close -text \"Close\" \\\n",\
"	    -command { \n",\
"		wm withdraw .memuse_dlg\n",\
"		set memuseplot 0\n",\
"	    }\n",\
"	pack $w.close\n",\
"	\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw $w\n",\
"	    wm geom $w +100+100\n",\
"	    wm deiconify $w\n",\
"	    wm title $w \"Memory Usage\"\n",\
"	    focus $w\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc STLinfodialog { show } {\n",\
"\n",\
"    set w .STLinfo_dlg\n",\
"    \n",\
"    if {[winfo exists .STLinfo_dlg] == 1} {\n",\
"\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw .STLinfo_dlg\n",\
"	    wm deiconify $w\n",\
"	    focus $w \n",\
"	} {\n",\
"	    wm withdraw $w\n",\
"	}\n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	set c $w.c\n",\
"\n",\
"	canvas $c -relief raised -width 450 -height 300\n",\
"	pack $w.c -side top -fill x\n",\
"\n",\
"	set plotFont {Helvetica 18}\n",\
"	set smallFont {Helvetica 12}\n",\
"\n",\
"	$c create line 100 250 400 250 -width 2\n",\
"	$c create line 100 250 100 50 -width 2\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu\n",\
"		\n",\
"	button $w.close -text \"Close\" \\\n",\
"	    -command { \n",\
"		wm withdraw .STLinfo_dlg\n",\
"			    }\n",\
"	pack $w.close\n",\
"    \n",\
"	\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw $w\n",\
"	    wm geom $w +100+100\n",\
"	    wm deiconify $w\n",\
"	    wm title $w \"STL Geometry Info\"\n",\
"	    focus $w\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc editprimitivedialog2 { name } {\n",\
"\n",\
"    global w classname\n",\
"\n",\
"    set w .ep_dlg\n",\
"    toplevel .$w\n",\
"\n",\
"    Ng_GetPrimitiveData $name classname valuelist\n",\
"        \n",\
"    \n",\
"    label $w.lab1 -text \"Primitive Name:  $name\";\n",\
"    label $w.lab2 -text \"Primitive Class: $classname\";\n",\
"    pack $w.lab1 $w.lab2 -fill x -pady 1m -padx 5m \n",\
"    \n",\
"    frame $w.specific -relief groove\n",\
"\n",\
"    global spec\n",\
"    set spec(sphere) { cx cy cz rad }\n",\
"    set spec(cylinder) { ax ay az bx by bz rad }\n",\
"    set spec(plane) { px py pz nx ny nz }\n",\
"    set spec(cone) { ax ay az bx by bz ra rb }\n",\
"    set spec(brick) { p1x p1y p1z p2x p2y p2z p3x p3y p3z p4x p4y p4z } \n",\
"   \n",\
"    set cnt 0\n",\
"    foreach field $spec($classname) {\n",\
"\n",\
"	frame $w.specific.f$cnt \n",\
"	pack $w.specific.f$cnt -side top -anchor ne\n",\
"\n",\
"	label $w.specific.f$cnt.lab -text \"$field\"\n",\
"	entry $w.specific.f$cnt.ent -textvariable dataval($cnt) \\\n",\
"	    -width 6 -relief sunken\n",\
"	pack $w.specific.f$cnt.ent $w.specific.f$cnt.lab -side right\n",\
"	$w.specific.f$cnt.ent delete 0 end\n",\
"	$w.specific.f$cnt.ent insert 0 [lindex $valuelist $cnt]\n",\
"	set cnt [expr $cnt + 1]\n",\
"    }\n",\
"    pack $w.specific\n",\
"\n",\
"\n",\
"    button $w.cancel -text \"cancel\" -command {\n",\
"	destroy $w \n",\
"    }\n",\
"\n",\
"    button $w.ok -text \"ok\" -command {\n",\
"\n",\
"	set valuelist \"\"\n",\
"	set cnt 0\n",\
"	foreach field $spec($classname) {\n",\
"	    lappend valuelist $dataval($cnt)\n",\
"	    set cnt [expr $cnt + 1]\n",\
"	}\n",\
"	Ng_SetPrimitiveData $name $valuelist\n",\
"	destroy $w\n",\
"    }\n",\
"    pack  $w.cancel $w.ok -side left -expand yes\n",\
"\n",\
"    bind $w <Return> { $w.ok  invoke}\n",\
"    bind $w <Escape> { $w.cancel  invoke}\n",\
"    \n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"\n",\
"    focus $w.specific.f0.ent\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc editprimitivedialog { } {\n",\
"    global w\n",\
"\n",\
"    set w .ep_dlg\n",\
"    toplevel $w\n",\
"\n",\
"    frame $w.frame -borderwidth 5m\n",\
"    pack $w.frame -side top -expand yes -fill y\n",\
"\n",\
"    listbox $w.frame.list -yscroll \"$w.frame.scroll set\" -setgrid 1 -height 12\n",\
"    scrollbar $w.frame.scroll -command \"$w.frame.list yview\"\n",\
"    pack $w.frame.scroll -side right -fill y\n",\
"    pack $w.frame.list -side left -expand 1 -fill both\n",\
"    \n",\
"\n",\
"    Ng_GetPrimitiveList primlist\n",\
"    foreach el $primlist {\n",\
"	$w.frame.list insert end $el }\n",\
"\n",\
"    button $w.cancel -text \"cancel\" -command { destroy $w }\n",\
"    button $w.ok -text \"ok\" -command {\n",\
"	set name [.ep_dlg.frame.list get active]\n",\
"	puts \"name=($name)\"\n",\
"	destroy $w\n",\
"	if { $name != \"\" } { editprimitivedialog2 $name }\n",\
"    }\n",\
"    \n",\
"    bind $w <Escape> { $w.cancel invoke }\n",\
"    bind $w <Return> { $w.ok invoke }\n",\
"    \n",\
"\n",\
"    pack  $w.cancel $w.ok -side left -expand yes\n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"\n",\
"    focus $w.frame.list\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc newprimitivedialog { } {\n",\
"\n",\
"    global w name\n",\
"\n",\
"    set w .ap_dlg\n",\
"    \n",\
"    toplevel $w\n",\
"\n",\
"    set name \"\"\n",\
"    frame $w.f1\n",\
"    pack $w.f1 -pady 2m\n",\
"    label $w.f1.lab -text \"Primitive Name: \";\n",\
"    entry $w.f1.ent -width 5 -relief sunken \\\n",\
"	-textvariable name\n",\
"    pack $w.f1.lab $w.f1.ent -side left\n",\
"    \n",\
"    frame $w.frame -borderwidth .5c\n",\
"    pack $w.frame -side top -expand yes -fill y\n",\
"\n",\
"    listbox $w.frame.list -yscroll \"$w.frame.scroll set\" -setgrid 1 -height 8 \n",\
"    scrollbar $w.frame.scroll -command \"$w.frame.list yview\"\n",\
"    pack $w.frame.scroll -side right -fill y\n",\
"    pack $w.frame.list -side left -expand 1 -fill both\n",\
"    \n",\
"    $w.frame.list insert 0 sphere cylinder plane cone brick\n",\
"    $w.frame.list activate 0\n",\
"    \n",\
"    button $w.ok -text \"ok\" -command {\n",\
"	Ng_CreatePrimitive [$w.frame.list get active]  $name\n",\
"	destroy $w\n",\
"	editprimitivedialog2 $name\n",\
"    }\n",\
"\n",\
"    button $w.cancel -text \"cancel\" -command {\n",\
"	destroy $w\n",\
"    }\n",\
"    \n",\
"    pack  $w.cancel $w.ok -side left -expand yes -pady 2m\n",\
"\n",\
"\n",\
"    bind $w <Escape> { $w.cancel invoke }\n",\
"    bind $w <Return> { $w.ok invoke }\n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"\n",\
"    focus $w.f1.ent\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc newsoliddialog { } {\n",\
"\n",\
"    global w name val sollist\n",\
"\n",\
"    set w .ns_dlg\n",\
"    toplevel $w\n",\
"\n",\
"    set name \"\"\n",\
"    frame $w.f1\n",\
"    label $w.f1.lab -text \"Solid Name: \";\n",\
"    entry $w.f1.ent -width 5 -relief sunken \\\n",\
"	-textvariable name\n",\
"    $w.f1.ent delete 0 end\n",\
"    button $w.f1.getsel -text \"Get Selected\" -command { \n",\
"	$w.f1.ent delete 0 end\n",\
"	$w.f1.ent insert 0 [$w.f3.list get active]\n",\
"	$w.bu.get invoke\n",\
"    }\n",\
"    pack $w.f1.getsel -side bottom\n",\
"    pack $w.f1.ent $w.f1.lab -side right\n",\
"\n",\
"\n",\
"    frame $w.f3 -borderwidth .5c\n",\
"    listbox $w.f3.list -yscroll \"$w.f3.scroll set\" -setgrid 1 -height 12\n",\
"    scrollbar $w.f3.scroll -command \"$w.f3.list yview\"\n",\
"    pack $w.f3.scroll -side right -fill y\n",\
"    pack $w.f3.list -side left -expand 1 -fill both\n",\
"    \n",\
"    Ng_GetSolidList sollist\n",\
"    foreach el $sollist {\n",\
"	$w.f3.list insert end $el }\n",\
"\n",\
"    frame $w.f2\n",\
"    label $w.f2.lab -text \"Solid Description: \";\n",\
"    pack $w.f2.lab\n",\
"\n",\
"\n",\
"    entry $w.f2.ent -width 100 -relief sunken \\\n",\
"	-textvariable val  -xscrollcommand \"$w.f2.scr set\"\n",\
"    scrollbar $w.f2.scr -relief sunken -orient horiz -command \\\n",\
"	\"$w.f2.ent xview\"\n",\
"    $w.f2.ent delete 0 end\n",\
"    pack $w.f2.ent $w.f2.scr -fill x\n",\
"\n",\
"\n",\
"\n",\
"    frame $w.bu\n",\
"    button $w.bu.close -text \"close\" -command {\n",\
"	destroy $w\n",\
"    }\n",\
"\n",\
"    button $w.bu.get -text \"get data\" -command {\n",\
"	Ng_GetSolidData $name val\n",\
"    }\n",\
"\n",\
"    button $w.bu.set -text \"set data\" -command {\n",\
"	Ng_SetSolidData $name $val\n",\
"    }\n",\
"\n",\
"    pack $w.bu.get $w.bu.set $w.bu.close -side left \n",\
"\n",\
"\n",\
"    pack $w.bu -pady 5 -side bottom                     ;    pack $w.f2 -pady 5 -side bottom                     ;    pack $w.f1 -pady 5 -side left                       ;    pack $w.f3 -side left -expand yes -fill y           ;\n",\
"\n",\
"\n",\
"    bind $w <Escape> { $w.bu.close invoke }\n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"\n",\
"    focus $w\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc toplevelproperties { w solname surfname } {\n",\
"\n",\
"    global properties\n",\
"\n",\
"    Ng_TopLevel getprop $solname $surfname properties\n",\
"\n",\
"\n",\
"    set w .tlprop_dlg\n",\
"\n",\
"    if {[winfo exists $w] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"    \n",\
"	label $w.lab1 -text \"Red\"\n",\
"	scale $w.scale1 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(red)\n",\
"	\n",\
"	label $w.lab2 -text \"Green\"\n",\
"	scale $w.scale2 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	-command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(green)\n",\
"	\n",\
"	label $w.lab3 -text \"Blue\"\n",\
"	scale $w.scale3 -orient horizontal -length 300 -from 0 -to 1 \\\n",\
"	    -resolution 0.01  -tickinterval 0.2 \\\n",\
"	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } -variable  properties(blue)\n",\
"\n",\
"	\n",\
"	pack $w.lab1 $w.scale1 $w.lab2 $w.scale2 $w.lab3 $w.scale3\n",\
"\n",\
"	checkbutton $w.cb4 -text \"Visible\" \\\n",\
"	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } \\\n",\
"	    -variable properties(visible)\n",\
"	\n",\
"	checkbutton $w.cb5 -text \"Transparent\" \\\n",\
"	    -command { Ng_TopLevel setprop $solname $surfname properties; redraw } \\\n",\
"	    -variable properties(transp)\n",\
"	\n",\
"	\n",\
"	pack  $w.cb4 $w.cb5\n",\
"	\n",\
"	\n",\
"	frame $w.bu\n",\
"	pack $w.bu -fill x\n",\
"	button $w.bu.ok -text \"Ok\" -command \"destroy .tlprop_dlg\"\n",\
"	pack $w.bu.ok  -expand yes\n",\
"    \n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	focus $w\n",\
"    }\n",\
"    wm title $w \"Properties $solname $surfname\"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc topleveldialog { } {\n",\
"\n",\
"    global w name val sollist\n",\
"\n",\
"    set w .tl_dlg\n",\
"    toplevel $w\n",\
"\n",\
"\n",\
"\n",\
"    frame $w.sol -borderwidth .5c\n",\
"    listbox $w.sol.list -yscroll \"$w.sol.scroll set\" -setgrid 1 -height 12\n",\
"    scrollbar $w.sol.scroll -command \"$w.sol.list yview\"\n",\
"    pack $w.sol.scroll -side right -fill y\n",\
"    pack $w.sol.list -side left -expand 1 -fill both\n",\
"    \n",\
"    Ng_GetSolidList sollist\n",\
"    foreach el $sollist {\n",\
"	$w.sol.list insert end $el }\n",\
"    Ng_GetPrimitiveList sollist\n",\
"    foreach el $sollist {\n",\
"	$w.sol.list insert end $el }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    frame $w.sul -borderwidth .5c\n",\
"    listbox $w.sul.list -yscroll \"$w.sul.scroll set\" -setgrid 1 -height 12\n",\
"    scrollbar $w.sul.scroll -command \"$w.sul.list yview\"\n",\
"    pack $w.sul.scroll -side right -fill y\n",\
"    pack $w.sul.list -side left -expand 1 -fill both\n",\
"    \n",\
"    Ng_GetSurfaceList sollist\n",\
"    foreach el $sollist {\n",\
"	$w.sul.list insert end $el }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    frame $w.topl -borderwidth .5c\n",\
"    listbox $w.topl.list -yscroll \"$w.topl.scroll set\" -setgrid 1 -height 12 \\\n",\
"	-command { puts hi }\n",\
"    scrollbar $w.topl.scroll -command \"$w.topl.list yview\"\n",\
"    pack $w.topl.scroll -side right -fill y\n",\
"    pack $w.topl.list -side left -expand 1 -fill both\n",\
"    \n",\
"    Ng_TopLevel getlist sollist\n",\
"    puts $sollist\n",\
"    foreach el $sollist {\n",\
"	set hel \"[ lindex $el 0 ]\"\n",\
"	if { [ llength $el ] == 2 } {\n",\
"	    set hel \"[ lindex $el 1 ] on [ lindex $el 0 ]\"\n",\
"	}\n",\
"	$w.topl.list insert end $hel \n",\
"    }\n",\
"\n",\
"\n",\
"    frame $w.bu\n",\
"\n",\
"    button $w.bu.close -text \"close\" -command {\n",\
"	destroy $w\n",\
"    }\n",\
"    button $w.bu.addsol -text \"Add Solid\" -command {\n",\
"	set solname [$w.sol.list get active]\n",\
"	Ng_TopLevel set $solname \"\"\n",\
"	Ng_ParseGeometry\n",\
"	$w.topl.list insert end $solname\n",\
"    }\n",\
"\n",\
"    button $w.bu.addsurf -text \"Add Surface\" -command {\n",\
"	set solname [$w.sol.list get active]\n",\
"	set surfname [$w.sul.list get active]\n",\
"	Ng_TopLevel set $solname $surfname\n",\
"	Ng_ParseGeometry\n",\
"	puts \"$solname on $surfname\"\n",\
"	$w.topl.list insert end \"$surfname on $solname\"\n",\
"    }\n",\
"\n",\
"    button $w.bu.remsol -text \"Remove\" -command {\n",\
"	set solname [$w.topl.list get active]\n",\
"	set surfname \"\"\n",\
"	if { [llength $solname] == 3 } {\n",\
"	    set surfname [lindex $solname 0]\n",\
"	    set solname [lindex $solname 2]\n",\
"	}\n",\
"	Ng_TopLevel remove $solname $surfname\n",\
"	Ng_ParseGeometry\n",\
"	$w.topl.list delete active\n",\
"    }\n",\
"\n",\
"    button $w.bu.prop -text \"Properties\" -command {\n",\
"	set solname [$w.topl.list get active]\n",\
"	set surfname \"\"\n",\
"	if { [llength $solname] == 3 } {\n",\
"	    set surfname [lindex $solname 0]\n",\
"	    set solname [lindex $solname 2]\n",\
"	}\n",\
"	toplevelproperties tlp $solname $surfname\n",\
"    }\n",\
"\n",\
"\n",\
"    \n",\
"\n",\
"    pack  $w.bu.close $w.bu.addsol $w.bu.addsurf $w.bu.remsol $w.bu.prop  -side left \n",\
"\n",\
"\n",\
"    pack $w.bu -side bottom\n",\
"    pack $w.sol -side left -expand yes -fill y           ;    pack $w.sul -side left -expand yes -fill y           ;    pack $w.topl -side left -expand yes -fill y           ;\n",\
"\n",\
"    bind $w <Escape> { $w.bu.close invoke }\n",\
"\n",\
"    wm withdraw $w\n",\
"    wm geom $w +100+100\n",\
"    wm deiconify $w\n",\
"\n",\
"    focus $w\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc topleveldialog2 { } {\n",\
"    set w .tl2_dlg\n",\
"    \n",\
"    if {[winfo exists .tl2_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	global name val sollist\n",\
"\n",\
"	frame $w.topl -borderwidth .5c\n",\
"	listbox $w.topl.list -yscroll \"$w.topl.scroll set\" -setgrid 1 -height 12\n",\
"	scrollbar $w.topl.scroll -command \"$w.topl.list yview\"\n",\
"	pack $w.topl.scroll -side right -fill y\n",\
"	pack $w.topl.list -side left -expand 1 -fill both\n",\
"	\n",\
"	Ng_TopLevel getlist sollist\n",\
"	puts $sollist\n",\
"	set i 1\n",\
"	foreach el $sollist {\n",\
"	    set hel \"$i: [ lindex $el 0 ]\"\n",\
"	    if { [ llength $el ] == 2 } {\n",\
"		set hel \"$i: [ lindex $el 1 ] on [ lindex $el 0 ]\"\n",\
"	    }\n",\
"	    incr i\n",\
"	    $w.topl.list insert end $hel }\n",\
"	\n",\
"	\n",\
"	frame $w.bu\n",\
"	\n",\
"	button $w.bu.close -text \"close\" -command {\n",\
"	    destroy .tl2_dlg\n",\
"	}\n",\
"	\n",\
"\n",\
"	button $w.bu.prop -text \"Properties\" -command {\n",\
"	    set solname [.tl2_dlg.topl.list get active]\n",\
"	    set surfname \"\"\n",\
"	    if { [llength $solname] == 2 } {\n",\
"		set solname [lindex $solname 1]\n",\
"	    }\n",\
"	    if { [llength $solname] == 4 } {\n",\
"		set surfname [lindex $solname 1]\n",\
"		set solname [lindex $solname 3]\n",\
"	    }\n",\
"	    toplevelproperties tlp $solname $surfname\n",\
"	}\n",\
"	\n",\
"	pack $w.bu.close $w.bu.prop  -side left \n",\
"	pack $w.bu -side bottom\n",\
"	pack $w.topl -side left -expand yes -fill y           ;		\n",\
"	bind .tl2_dlg.topl.list <Double-1> {\n",\
"	    set solname [.tl2_dlg.topl.list get  @%x,%y]\n",\
"	    set surfname \"\"\n",\
"	    if { [llength $solname] == 2 } {\n",\
"		set solname [lindex $solname 1]\n",\
"	    }\n",\
"	    if { [llength $solname] == 4 } {\n",\
"		set surfname [lindex $solname 1]\n",\
"		set solname [lindex $solname 3]\n",\
"	    }\n",\
"	    toplevelproperties tlp $solname $surfname\n",\
"	}\n",\
"	\n",\
"	bind .tl2_dlg <Escape> { .tl2_dlg.bu.close invoke }\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Top-Level Options\"	\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc logwindow { } {\n",\
"    set w .logwindow\n",\
"    \n",\
"    if {[winfo exists .logwindow] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	text $w.edit -yscroll \"$w.scrolly set\" -setgrid 1 -height 12\n",\
"	scrollbar $w.scrolly -command \"$w.edit yview\"	\n",\
"	pack $w.edit -side left -fill both -expand 1\n",\
"	pack $w.scrolly -side left -fill both -expand 0\n",\
"\n",\
"	.logwindow.edit insert end \"Netgen Log Window\\n\"\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Netgen Log\"	\n",\
"	focus $w\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"set entities [ ]\n",\
"\n",\
"\n",\
"proc occdialogbuildtree {} {\n",\
"    global entities\n",\
"\n",\
"    set w .occ_dlg\n",\
"    set hlist [$w.mtre subwidget hlist]\n",\
"\n",\
"    set entities [Ng_GetOCCData getentities]\n",\
"    set nrentities [expr [llength $entities]]\n",\
"\n",\
"\n",\
"    if {$nrentities != 0} {\n",\
"\n",\
"	$hlist add Topology -itemtype text -text \"Topology\"\n",\
"	\n",\
"	$hlist add Topology/CompSolids   -itemtype text -text \"Composite Solids\" -data \"Composite Solids\"\n",\
"	$hlist add Topology/FreeSolids   -itemtype text -text \"Free Solids\" -data \"Free Solids\"\n",\
"	$hlist add Topology/FreeShells   -itemtype text -text \"Free Shells\" -data \"Free Shells\"\n",\
"	$hlist add Topology/FreeFaces    -itemtype text -text \"Free Faces\" -data \"Free Faces\"\n",\
"	$hlist add Topology/FreeWires    -itemtype text -text \"Free Wires\" -data \"Free Wires\"\n",\
"	$hlist add Topology/FreeEdges    -itemtype text -text \"Free Edges\" -data \"Free Edges\"\n",\
"	$hlist add Topology/FreeVertices -itemtype text -text \"Free Vertices\" -data \"Free Vertices\"\n",\
"\n",\
"	\n",\
"	set i [expr 0]\n",\
"	while {$i < $nrentities} {\n",\
"	    set entity [lindex $entities [expr $i]]\n",\
"	    incr i 1\n",\
"	    set entityname [lindex $entities [expr $i]]\n",\
"	    $hlist add Topology/$entity -text $entityname -data $entityname\n",\
"	    incr i 1\n",\
"	    $w.mtre close Topology/$entity\n",\
"	}\n",\
"	\n",\
"	$w.mtre autosetmode\n",\
"	\n",\
"	$w.mtre open Topology\n",\
"	$w.mtre close Topology/CompSolids\n",\
"	$w.mtre close Topology/FreeSolids\n",\
"	$w.mtre close Topology/FreeShells\n",\
"	$w.mtre close Topology/FreeFaces\n",\
"	$w.mtre close Topology/FreeWires\n",\
"	$w.mtre close Topology/FreeEdges\n",\
"	$w.mtre close Topology/FreeVertices\n",\
"	\n",\
"	set i [expr 0]\n",\
"	while {$i < $nrentities} {\n",\
"	    set entity [lindex $entities [expr $i]]\n",\
"	    $w.mtre close Topology/$entity\n",\
"	    incr i 2\n",\
"	}\n",\
"	\n",\
"	set faces [Ng_OCCCommand getunmeshedfaceinfo]    \n",\
"	set nrfaces [expr [llength $faces]]\n",\
"	if {$nrfaces >= 2} {\n",\
"	    $hlist add ErrorFaces -itemtype text -text \"Faces with surface meshing error\"\n",\
"	    $w.mtre open ErrorFaces\n",\
"	    set i [expr 0]\n",\
"	    while {$i < $nrfaces} {\n",\
"		set entity [lindex $faces [expr $i]]\n",\
"		incr i 1\n",\
"		set entityname [lindex $faces [expr $i]]\n",\
"		$hlist add ErrorFaces/$entity -text $entityname -data $entityname\n",\
"		incr i 1\n",\
"	    }\n",\
"	}\n",\
"	\n",\
"\n",\
"	set faces [Ng_OCCCommand getnotdrawablefaces]    \n",\
"	set nrfaces [expr [llength $faces]]\n",\
"	if {$nrfaces >= 2} {\n",\
"	    $hlist add NotDrawableFaces -itemtype text -text \"Faces impossible to visualize\"\n",\
"	    $w.mtre open NotDrawableFaces\n",\
"	    set i [expr 0]\n",\
"	    while {$i < $nrfaces} {\n",\
"		set entity [lindex $faces [expr $i]]\n",\
"		incr i 1\n",\
"		set entityname [lindex $faces [expr $i]]\n",\
"		$hlist add NotDrawableFaces/$entity -text $entityname -data $entityname\n",\
"		incr i 1\n",\
"	    }\n",\
"	}\n",\
"\n",\
"\n",\
"	$w.mtre autosetmode\n",\
"\n",\
"	puts \"done\"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"proc rebuildoccdialog {} {\n",\
"    if {[winfo exists .occ_dlg] == 1} {\n",\
"	[.occ_dlg.mtre subwidget hlist] delete all\n",\
"	occdialogbuildtree \n",\
"    }\n",\
"}\n",\
"\n",\
"proc checkoccloaded { } {\n",\
"    set isoccgeometryloaded [Ng_OCCCommand isoccgeometryloaded]\n",\
"    if {$isoccgeometryloaded == 0} {\n",\
"	puts \"no IGES/STEP geometry loaded\"\n",\
"	destroy .occ_dlg\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc selectentity { entityname } {\n",\
"    global entities\n",\
"    set nrentities [expr [llength $entities]]\n",\
"    set i [expr 0]\n",\
"    while {$i < $nrentities} {\n",\
"	set entitylength []\n",\
"			\n",\
"	set entity2 [lindex $entities [expr $i]]\n",\
"	incr i 1\n",\
"	set entityname2 [lindex $entities [expr $i]]\n",\
"	incr i 1\n",\
"	set entityname2 [string range $entityname2 0 [expr [string length $entityname]-1]]\n",\
"	\n",\
"	if {$entityname == $entityname2} {\n",\
"	    set hlist [.occ_dlg.mtre subwidget hlist]\n",\
"	    .occ_dlg.mtre open Topology\n",\
"	    set slashpos [string last \"/\" $entity2]\n",\
"	    set entity3 [string range $entity2 0 [expr $slashpos-1]]\n",\
"	    while {$slashpos != -1} {\n",\
"		.occ_dlg.mtre open Topology/$entity3\n",\
"		\n",\
"		set slashpos [string last \"/\" $entity3]\n",\
"		set entity3 [string range $entity3 0 [expr $slashpos-1]]\n",\
"	    }\n",\
"	    $hlist selection clear\n",\
"	    $hlist see Topology/$entity2\n",\
"	    $hlist selection set Topology/$entity2\n",\
"	} \n",\
"    }	    \n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc occdialog { } {\n",\
"	\n",\
"    uplevel 1 {\n",\
"\n",\
"    global entities\n",\
"    set selectvisual geometry\n",\
"    Ng_SetVisParameters\n",\
"    redraw\n",\
"\n",\
"    set w .occ_dlg\n",\
"    \n",\
"    if {[winfo exists .occ_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {	\n",\
"	toplevel $w\n",\
"\n",\
"	tixTree $w.mtre -options { separator \"/\" }\n",\
"	pack $w.mtre -fill both -expand yes\n",\
"\n",\
"	occdialogbuildtree\n",\
"\n",\
"	set hlist [$w.mtre subwidget hlist]\n",\
"\n",\
"\n",\
"	set solname {\"\"}\n",\
"\n",\
"	\n",\
"	bind $hlist <Double-1> {\n",\
"	    set oldsolname {$solname}\n",\
"	    set solname [[.occ_dlg.mtre subwidget hlist] info selection]\n",\
"	    if {$solname != \"\" && $oldsolname != $solname } {\n",\
"		set seppos [string first \"/\" $solname]\n",\
"		set rootname [string range $solname 0 [expr $seppos-1]]\n",\
"\n",\
"		set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]\n",\
"		set spacepos [string first \" \" $entityname]\n",\
"		set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"		set spacepos2 [string first \" \" $helpstring]\n",\
"		set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"		if {$rootname == \"Topology\"} {\n",\
"		    Ng_OCCCommand highlightentity $entitytype $entitynumber\n",\
"		    set selectvisual geometry\n",\
"		    redraw\n",\
"		} {\n",\
"		    set brackpos [string first \" (\" $entityname]\n",\
"		    if {$brackpos != -1} {\n",\
"			set entityname [string range $entityname 0 $brackpos]\n",\
"		    }\n",\
"\n",\
"		    selectentity $entityname\n",\
"		}\n",\
"	    }\n",\
"	}\n",\
"	\n",\
"	button $w.cl -text \"Close\" -command {\n",\
"	    destroy .occ_dlg\n",\
"	}\n",\
"	\n",\
"	button $w.show -text \"Show\" -command {\n",\
"	    set solname [[.occ_dlg.mtre subwidget hlist] info selection]\n",\
"	    set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]\n",\
"	    set spacepos [string first \" \" $entityname]\n",\
"	    set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"	    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"	    set spacepos2 [string first \" \" $helpstring]\n",\
"	    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"	    Ng_OCCCommand show $entitytype $entitynumber\n",\
"	    set selectvisual geometry\n",\
"	    redraw\n",\
"	}\n",\
"	button $w.hide -text \"Hide\" -command {\n",\
"	    set solname [[.occ_dlg.mtre subwidget hlist] info selection]\n",\
"	    set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]\n",\
"	    set spacepos [string first \" \" $entityname]\n",\
"	    set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"	    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"	    set spacepos2 [string first \" \" $helpstring]\n",\
"	    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"	    Ng_OCCCommand hide $entitytype $entitynumber\n",\
"	    set selectvisual geometry\n",\
"	    redraw\n",\
"	}\n",\
"\n",\
"	button $w.swaporientation -text \"Swap orientation\" -command {\n",\
"	    set solname [[.occ_dlg.mtre subwidget hlist] info selection]\n",\
"	    set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]\n",\
"	    set spacepos [string first \" \" $entityname]\n",\
"	    set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"	    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"	    set spacepos2 [string first \" \" $helpstring]\n",\
"	    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"	    Ng_OCCCommand swaporientation $entitytype $entitynumber\n",\
"	    set selectvisual geometry\n",\
"	    redraw\n",\
"\n",\
"	    [.occ_dlg.mtre subwidget hlist] delete all\n",\
"	    occdialogbuildtree	\n",\
"	}\n",\
"\n",\
"	button $w.marksingular -text \"Mark/Unmark as singular\" -command {\n",\
"	    set solname [[.occ_dlg.mtre subwidget hlist] info selection]\n",\
"	    set entityname [[.occ_dlg.mtre subwidget hlist] info data $solname]\n",\
"	    set spacepos [string first \" \" $entityname]\n",\
"	    if { $spacepos != 0 } {\n",\
"		set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"		set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"		set spacepos2 [string first \" \" $helpstring]\n",\
"		if { $spacepos2 != 0 } {\n",\
"		    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"		    \n",\
"		    global ismarkedsingular\n",\
"		    Ng_OCCCommand marksingular $entitytype $entitynumber\n",\
"		    \n",\
"		    set hlist [$w.mtre subwidget hlist]\n",\
"		    \n",\
"		    		    		    set style1 [tixDisplayStyle imagetext -foreground black -background white -selectforeground white -selectbackground blue]\n",\
"		    set style2 [tixDisplayStyle imagetext -foreground red -background white -selectforeground red -selectbackground blue]\n",\
"		    \n",\
"		    if { $ismarkedsingular == 0 } {\n",\
"			$hlist entryconfigure $solname -style $style1\n",\
"		    } {\n",\
"			$hlist entryconfigure $solname -style $style2\n",\
"		    }\n",\
"\n",\
"		}\n",\
"	    }\n",\
"\n",\
"\n",\
"	}\n",\
"\n",\
"\n",\
"	checkbutton $w.zoomtohighlightedentity -text \"Zoom to highlighted entity\" \\\n",\
"	    -variable occoptions.zoomtohighlightedentity \\\n",\
"	    -command {\n",\
"		Ng_SetOCCVisParameters\n",\
"		if { ${occoptions.zoomtohighlightedentity} == 1} {\n",\
"		    set selectvisual geometry\n",\
"		    Ng_OCCCommand redrawstatus 1\n",\
"		    redraw\n",\
"		} {\n",\
"		    Ng_OCCCommand redrawstatus 0\n",\
"		}\n",\
"	    }\n",\
"\n",\
"\n",\
"\n",\
"	frame $w.healing -relief groove -borderwidth 3\n",\
"\n",\
"	button $w.healing.checkentities -text \"Analyze geometry\" -command {\n",\
"	    set irregent [Ng_OCCCommand findsmallentities]\n",\
"\n",\
"	    set w .occ_dlg\n",\
"	    set hlist [$w.mtre subwidget hlist]\n",\
"	    \n",\
"	    $hlist add ProblematicEntities -text \"Problematic Entities\"\n",\
"	    $hlist delete offsprings ProblematicEntities\n",\
"\n",\
"	    set nritems [expr [llength $irregent]]\n",\
"	    set i [expr 0]\n",\
"	    while {$i < $nritems} {\n",\
"		set entity [lindex $irregent [expr $i]]\n",\
"		incr i 1\n",\
"		set entityname [lindex $irregent [expr $i]]\n",\
"		$hlist add ProblematicEntities/$entity -text $entityname -data $entityname\n",\
"		incr i 1\n",\
"	    }\n",\
"	    $w.mtre open ProblematicEntities\n",\
"	    $w.mtre autosetmode\n",\
"	}\n",\
"\n",\
"	tixControl $w.healing.tolerance -label \"Healing tolerance: \" -integer false \\\n",\
"	    -variable occoptions.tolerance -min 1e-9 -max 1e6 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	checkbutton $w.healing.fixsmalledges -text \"Fix small edges\" \\\n",\
"	    -variable occoptions.fixsmalledges\n",\
"	\n",\
"	checkbutton $w.healing.fixspotstripfaces -text \"Fix spot/strip faces\" \\\n",\
"	    -variable occoptions.fixspotstripfaces\n",\
"	\n",\
"	checkbutton $w.healing.sewfaces -text \"Sew faces\" \\\n",\
"	    -variable occoptions.sewfaces\n",\
"	\n",\
"	checkbutton $w.healing.makesolids -text \"Make solids\" \\\n",\
"	    -variable occoptions.makesolids\n",\
"	\n",\
"	checkbutton $w.healing.splitpartitions -text \"Split partitions\" \\\n",\
"	    -variable occoptions.splitpartitions\n",\
"	\n",\
"	button $w.healing.heal -text \"Heal geometry\" -command { \n",\
"	    .occ_dlg.healing.tolerance invoke\n",\
"	    Ng_OCCCommand shapehealing\n",\
"	    redraw \n",\
"	    [.occ_dlg.mtre subwidget hlist] delete all\n",\
"	    occdialogbuildtree\n",\
"	}\n",\
"\n",\
"	pack $w.healing.checkentities\n",\
"\n",\
"	pack $w.healing.tolerance $w.healing.fixsmalledges \\\n",\
"	    $w.healing.fixspotstripfaces $w.healing.sewfaces \\\n",\
"	    $w.healing.makesolids $w.healing.splitpartitions -anchor w\n",\
"\n",\
"	pack $w.healing.heal	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	pack $w.show $w.hide\n",\
"\n",\
"	pack $w.zoomtohighlightedentity -anchor w\n",\
"	pack $w.swaporientation\n",\
"	pack $w.marksingular\n",\
"	pack $w.healing -fill x\n",\
"	pack $w.cl\n",\
" \n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"IGES/STEP Topology Explorer/Doctor\"\n",\
"	focus .occ_dlg\n",\
"    }\n",\
"}\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc printtable { tablevar } {\n",\
"    set w newtcltable\n",\
"    while {[winfo exists .$w] == 1} {set w 1$w}\n",\
"    set w .$w\n",\
"    toplevel $w\n",\
"    for {set i 0} {$i < [lindex $tablevar 2]} { incr i } {\n",\
"	frame $w.col$i\n",\
"	for {set j 0} {$j < [lindex $tablevar 1]} { incr j } {\n",\
"	    frame $w.col$i.row$j\n",\
"	    message $w.col$i.row$j.txt -aspect 10000000 -text [lindex $tablevar [expr 3+[lindex $tablevar 2]*$j+$i]]\n",\
"	    pack $w.col$i.row$j.txt\n",\
"	    pack $w.col$i.row$j -side top\n",\
"	}\n",\
"	pack $w.col$i -side left\n",\
"    }\n",\
"    wm withdraw $w\n",\
"    wm geom $w +200+100; wm deiconify $w\n",\
"    wm title $w [lindex $tablevar 0]\n",\
"    focus $w\n",\
"}\n",\
"\n",\
"\n",\
"set latestwarning 0\n",\
"\n",\
"\n",\
"proc printwarning { textvar } {\n",\
"    global latestwarning\n",\
"    set latestwarning $textvar\n",\
"    set w warning\n",\
"    while {[winfo exists .$w] == 1} {set w 1$w}\n",\
"    set w .$w\n",\
"    toplevel $w\n",\
"    message $w.mes -aspect 2000 -text \"WARNING:\\n$textvar\"\n",\
"    button $w.done -text \"Done\" -command \"destroy $w\"\n",\
"    pack $w.mes\n",\
"    pack $w.done\n",\
"    wm withdraw $w\n",\
"    wm deiconify $w\n",\
"    wm title $w \"Warning\"\n",\
"    focus $w\n",\
"}\n",\
"\n",\
"\n",\
"proc printlatestwarning { } { \n",\
"    global latestwarning \n",\
"    if {$latestwarning != 0} {printwarning $latestwarning}\n",\
"}\n",\
"\n",\
"\n",\
"proc paralleldialog { } {\n",\
"\n",\
"    set w .parallel_dlg\n",\
"    \n",\
"    if {[winfo exists .parallel_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	wm geometry $w =270x100\n",\
"\n",\
"	focus $w \n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"	wm geometry $w =270x100\n",\
"\n",\
"	set ww $w\n",\
"\n",\
"	button $ww.visallb -text \"View All\"   -width 20 -command\\\n",\
"	    { Ng_VisualizeAll; } \n",\
"	pack $ww.visallb  \n",\
"	\n",\
"	button $ww.visoneb -text \"View One\"   -width 20 -command \\\n",\
"	    { Ng_VisualizeOne; } \n",\
"	pack $ww.visoneb  \n",\
"	\n",\
"	button $ww.overlap -text \"overlap++\"    -width 20 -command \\\n",\
"	    { Ng_IncrOverlap; }\n",\
"	\n",\
"	pack $ww.overlap \n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Parallel Netgen\"\n",\
"	focus .parallel_dlg \n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc runtestdialog { } {\n",\
"    source $::ngdir/ngtcltk/ngshell.tcl\n",\
"    set w .runtest_dlg\n",\
"    \n",\
"    if {[winfo exists .runtest_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	frame $w.in2dframe \n",\
"	pack $w.in2dframe\n",\
"\n",\
"        set in2dlogfile \"\"\n",\
" 	tixLabelEntry $w.in2dframe.ent -label \"in2d log-file: console if empty\"  \\\n",\
" 	    -labelside top \\\n",\
" 	    -options {  \n",\
" 		entry.textVariable in2dlogfile\n",\
" 		entry.width 35\n",\
" 		label.width 25\n",\
" 		label.anchor w\n",\
" 	    }	\n",\
" 	button $w.in2dframe.btn -text \"Browse\" -command {\n",\
"	    set types { { \"Log file\"   {.log}	} }\n",\
"	    set in2dlogfile [tk_getOpenFile -filetypes $types -initialfile $in2dlogfile]\n",\
"	}\n",\
" 	button $w.in2dframe.test -text \"Test in2d meshing\" -command { ngtest in2d $in2dlogfile }\n",\
"\n",\
"        \n",\
" 	pack $w.in2dframe.test -side left -anchor s -padx 4 -pady 4\n",\
" 	pack $w.in2dframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4\n",\
" 	pack $w.in2dframe.btn -side left -anchor s -padx 4 -pady 4\n",\
"\n",\
"        \n",\
"	frame $w.geoframe \n",\
"	pack $w.geoframe\n",\
"\n",\
"        set geologfile \"\" \n",\
" 	tixLabelEntry $w.geoframe.ent -label \"geo log-file: console if empty\"  \\\n",\
" 	    -labelside top \\\n",\
" 	    -options {  \n",\
" 		entry.textVariable geologfile\n",\
" 		entry.width 35\n",\
" 		label.width 25\n",\
" 		label.anchor w\n",\
" 	    }	\n",\
" 	button $w.geoframe.btn -text \"Browse\" -command {\n",\
"	    set types { { \"Log file\"   {.log}	} }\n",\
"	    set geologfile [tk_getOpenFile -filetypes $types -initialfile $geologfile]\n",\
"	}\n",\
" 	button $w.geoframe.test -text \"Test geo meshing\" -command { ngtest geo $geologfile }\n",\
"\n",\
"        \n",\
" 	pack $w.geoframe.test -side left -anchor s -padx 4 -pady 4\n",\
" 	pack $w.geoframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4\n",\
" 	pack $w.geoframe.btn -side left -anchor s -padx 4 -pady 4\n",\
"\n",\
"	frame $w.stlframe \n",\
"	pack $w.stlframe\n",\
"\n",\
"        set stllogfile \"\"\n",\
" 	tixLabelEntry $w.stlframe.ent -label \"stl log-file: console if empty\"  \\\n",\
" 	    -labelside top \\\n",\
" 	    -options {  \n",\
" 		entry.textVariable stllogfile\n",\
" 		entry.width 35\n",\
" 		label.width 25\n",\
" 		label.anchor w\n",\
" 	    }	\n",\
" 	button $w.stlframe.btn -text \"Browse\" -command {\n",\
"	    set types { { \"Log file\"   {.log}	} }\n",\
"	    set stllogfile [tk_getOpenFile -filetypes $types -initialfile $stllogfile]\n",\
"	}\n",\
" 	button $w.stlframe.test -text \"Test stl meshing\" -command { ngtest stl $stllogfile }\n",\
"\n",\
"        \n",\
" 	pack $w.stlframe.test -side left -anchor s -padx 4 -pady 4\n",\
" 	pack $w.stlframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4\n",\
" 	pack $w.stlframe.btn -side left -anchor s -padx 4 -pady 4\n",\
"\n",\
"	frame $w.pdeframe \n",\
"	pack $w.pdeframe\n",\
"\n",\
"        set pdelogfile \"\"\n",\
" 	tixLabelEntry $w.pdeframe.ent -label \"pde log-file: console if empty\"  \\\n",\
" 	    -labelside top \\\n",\
" 	    -options {  \n",\
" 		entry.textVariable pdelogfile\n",\
" 		entry.width 35\n",\
" 		label.width 25\n",\
" 		label.anchor w\n",\
" 	    }	\n",\
" 	button $w.pdeframe.btn -text \"Browse\" -command {\n",\
"	    set types { { \"Log file\"   {.log}	} }\n",\
"	    set pdelogfile [tk_getOpenFile -filetypes $types -initialfile $pdelogfile]\n",\
"	}\n",\
" 	button $w.pdeframe.test -text \"Test ngsolve pde's\" -command { ngtest pde $pdelogfile }\n",\
"\n",\
"        \n",\
" 	pack $w.pdeframe.test -side left -anchor s -padx 4 -pady 4\n",\
" 	pack $w.pdeframe.ent -side left -expand yes -fill x -anchor s -padx 4 -pady 4\n",\
" 	pack $w.pdeframe.btn -side left -anchor s -padx 4 -pady 4\n",\
" \n",\
"	wm title $w \"Testing\"\n",\
"	focus .runtest_dlg \n",\
"    }\n",\
"}\n",\
"\n",\
"set oldmousex 0\n",\
"set oldmousey 0\n",\
"if {[catch {togl .ndraw -width 400 -height 300  -rgba true -double true -depth true -privatecmap false -stereo false -indirect false }] } {    \n",\
"    puts \"no OpenGL\" \n",\
"} {\n",\
"        pack .ndraw -expand true -fill both -padx 10 -pady 10\n",\
"        bind .ndraw <Button-1> {\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"    bind .ndraw <Button-2> {\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"    bind .ndraw <Button-3> {\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"    bind .ndraw <B1-Motion> {\n",\
"	Ng_MouseMove $oldmousex $oldmousey %x %y $drawmode\n",\
"	.ndraw render\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"\n",\
"    bind .ndraw <Double-1> {\n",\
"	Ng_MouseDblClick %x %y\n",\
"	.ndraw render\n",\
"	if { [winfo exists .bcprop_dlg] } { bcpropdialog }\n",\
"	if { [winfo exists .fieldlines_dlg] } { fieldlinesdialog }\n",\
"    }\n",\
"\n",\
"    bind .ndraw <B2-Motion> {\n",\
"	Ng_MouseMove $oldmousex $oldmousey %x %y move\n",\
"	.ndraw render\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"\n",\
"    bind .ndraw <B3-Motion> {\n",\
"	Ng_MouseMove $oldmousex $oldmousey %x %y zoom\n",\
"	.ndraw render\n",\
"	set oldmousex %x; set oldmousey %y;\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"proc popupcheckredraw { vari { x 0 } } {\n",\
"    upvar $vari varname\n",\
"    if { $varname == 1 } {\n",\
"	set varname 0\n",\
"    } {\n",\
"        	Ng_Vis_Set parameters\n",\
"	redraw\n",\
"    }\n",\
"}\n",\
"proc popupcheckredraw2 { vari boolvar { x 0 } } {\n",\
"    upvar $vari varname\n",\
"    if { $varname == 1 } {\n",\
"	set varname 0\n",\
"    } {\n",\
"	Ng_SetVisParameters\n",\
"	if { $boolvar == 1 } { redraw }\n",\
"	Ng_SetVisParameters\n",\
"    }\n",\
"}\n",\
"proc popupcheckredraw3 { vari { x 0 } } {\n",\
"    upvar $vari varname\n",\
"    if { $varname == 1 } {\n",\
"	set varname 0\n",\
"    } {\n",\
"	Ng_Vis_Set parameters\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc redraw { {x 0} } {\n",\
"    if {[winfo exists .ndraw]} { .ndraw render } \n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"bind . <Left> { Ng_MouseMove 0 0 -10 0 rotate; redraw }\n",\
"bind . <Right> { Ng_MouseMove 0 0 10 0 rotate; redraw }\n",\
"bind . <Up> { Ng_MouseMove 0 0 0 -10 rotate; redraw }\n",\
"bind . <Down> { Ng_MouseMove 0 0 0 10 rotate; redraw }\n",\
"bind . <Shift-Left> { Ng_MouseMove 0 0 -10 0 move; redraw }\n",\
"bind . <Shift-Right> { Ng_MouseMove 0 0 10 0 move; redraw }\n",\
"bind . <Shift-Up> { Ng_MouseMove 0 0 0 -10 move; redraw }\n",\
"bind . <Shift-Down> { Ng_MouseMove 0 0 0 10 move; redraw }\n",\
"bind . <Control-Up> { Ng_MouseMove 0 0 0 -10 zoom; redraw }\n",\
"bind . <Control-Down> { Ng_MouseMove 0 0 0 10 zoom; redraw }\n",\
"\n",\
"bind all <Button-4> \\\n",\
"   {event generate [focus -displayof %W] <MouseWheel> -delta  120}\n",\
"\n",\
" bind all <Button-5> \\\n",\
"   {event generate [focus -displayof %W] <MouseWheel> -delta -120}\n",\
"\n",\
"bind all <MouseWheel> { Ng_MouseMove 0 0 0 [expr {%D/-5}] zoom; redraw }\n",\
"\n",\
"proc print_commandline_help { } {\n",\
"    \n",\
"    puts \"Usage: ng { options }\"\n",\
"\n",\
"    puts \"-geofile=filename      Input geometry file (alternative:  ng filename)\"\n",\
"    puts \"-meshfile=filename     Output mesh file\"\n",\
"    puts \"-verycoarse, -coarse, -moderate, -fine, -veryfine\"\n",\
"    puts \"                       Automatic mesh-size selection\"\n",\
"    puts \"-meshsizefile=filename Load mesh-size file with local mesh sizes\"\n",\
"    puts \"-meshfiletype={\\\"Neutral Format\\\", ...}\"\n",\
"    puts \"                       Filetype of output file, default is netgen file\"\n",\
"    puts \"-batchmode             Run Netgen in batchmode\"\n",\
"    puts \"-inputmeshfile=filename\"\n",\
"    puts \"                       Input mesh file (batchmode only)\"\n",\
"    puts \"-mergefile=filename    Merge with mesh file (batchmode only)\"\n",\
"    puts \"-refinementfile=filename\"\n",\
"    puts \"                       Use refinementinfo from file (batchmode only)\"\n",\
"    puts \"-serversocket=\\#num    Start a Netgen server with port \\#num\"\n",\
"    puts \"-V                     Print additional information\"\n",\
"    puts \"-testout=filename      file for test output\"\n",\
"\n",\
"    if { [catch { NGS_GetData } ] == 0 } { \n",\
"	puts \"\\nNGSolve parameters:\"\n",\
"	puts \"-pdefile=filename      Load pde input file\"\n",\
"	puts \"-solve                 Solve pde once\"\n",\
"	puts \"-solve=n               Solve pde by n adaptive refinement steps\"\n",\
"	puts \"-recent                Load and solve most recently loaded pde\"\n",\
"    }\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc set_menu_help { entry helpmsg } {\n",\
"    global menuhelps\n",\
"    set menuhelps($entry) $helpmsg\n",\
"}\n",\
"\n",\
"proc show_menu_help { entry } {\n",\
"    global menuhelps\n",\
"\n",\
"\n",\
"    if {[catch {set helptext $menuhelps($entry)}]} {\n",\
"	set helptext \"no help available   \"\n",\
"    }    \n",\
"\n",\
"    .helpline configure -text $helptext\n",\
"    \n",\
"    if {[winfo exists .senshelp_dlg]==1} {\n",\
"	.senshelp_dlg.text delete 1.0 end\n",\
"	.senshelp_dlg.text insert end \"Menu item: $entry\\n\\n\"\n",\
"	.senshelp_dlg.text insert end $helptext\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"tixBalloon .balloon -statusbar .helpline\n",\
"\n",\
"proc set_control_help { control helpmsg } {\n",\
"    bind $control <Enter> \"show_control_help {$helpmsg}\"\n",\
"    bind $control <Leave> \"show_control_help {None}\"\n",\
"    .balloon bind  $control -balloonmsg $helpmsg -statusmsg $helpmsg\n",\
"}\n",\
"\n",\
"proc show_control_help { helpmsg } {\n",\
"    .helpline configure -text $helpmsg\n",\
"    if {[winfo exists .senshelp_dlg]==1} {\n",\
"	.senshelp_dlg.text delete 1.0 end\n",\
"	.senshelp_dlg.text insert end $helpmsg\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"proc sensitivehelpdialog { show } {\n",\
"\n",\
"    set w .senshelp_dlg\n",\
"    \n",\
"    if {[winfo exists .senshelp_dlg] == 1} {\n",\
"\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw .senshelp_dlg\n",\
"	    wm deiconify $w\n",\
"	    focus $w \n",\
"	} {\n",\
"	    wm withdraw $w\n",\
"	}\n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	global senshelptext\n",\
"\n",\
"	text $w.text -yscrollcommand \"$w.scroll set\" -setgrid true \\\n",\
"	    -width 40 -height 10  -wrap word\n",\
"	scrollbar $w.scroll -command \"$w.text yview\"\n",\
"	pack $w.scroll -side right -fill y\n",\
"	pack $w.text -expand yes -fill both\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu\n",\
"		\n",\
"	button $w.close -text \"Close\" \\\n",\
"	    -command { \n",\
"		wm withdraw .senshelp_dlg\n",\
"		set showsensitivehelp 0\n",\
"	    }\n",\
"	pack $w.close\n",\
"    \n",\
"	\n",\
"	if { $show == 1 } {\n",\
"	    wm withdraw $w\n",\
"	    wm geom $w +100+100\n",\
"	    wm deiconify $w\n",\
"	    wm title $w \"Help\"\n",\
"	    focus $w\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"set_menu_help \"File\"  \"In File menu you can load and store geometries, meshes etc.\" \n",\
"\n",\
"set_menu_help \"New Geometry\"  \"Deletes current geometry\"\n",\
"set_menu_help \"Load Geometry\"  \"Loads Geometry file in one of the formats STL (ASCII or binary), Constructive Solid Geometry (.geo) or 2D geometry. Please have a look into Netgen User's manuel for more details.\"\n",\
"set_menu_help \"Save Geometry\" \"Saves STL Geometry in in either ASCII or binary STL format.\"\n",\
"set_menu_help \"Load Mesh\" \"Loads surface and volume mesh in Netgen internal format.\"\n",\
"set_menu_help \"Save Mesh\" \"Saves surface and volume mesh in Netgen internal format.\"\n",\
"set_menu_help \"Write EPS File\" \"Dumps OpenGL rendering to EPS File.\"\n",\
"set_menu_help \"Save Options\" \"Saves current options in file \\\"ng.opt\\\". These options will be loaded again when starting ng in the same directory.\"\n",\
"set_menu_help \"Export Mesh\" \"Exports mesh in format defined by Export Filetype.\"\n",\
"set_menu_help \"Export Filetype\" \"Selects file format for exporting mesh. Please have a look into the Netgen User's manual for more information.\"\n",\
"set_menu_help \"Import Mesh\" \"Imports surface or volume mesh in exchange format.\"\n",\
"set_menu_help \"Quit\" \"Quits Netgen\"\n",\
"\n",\
"set_menu_help \"Geometry\" \"Preparing geometries, visualiztion of geometries.\"\n",\
"set_menu_help \"Scan CSG Geometry\" \"Generates surface triangulation for rendering\"\n",\
"set_menu_help \"CSG Options\" \"Sets Options for CSG visualization (bounding box, detail size, number of facets).\"\n",\
"set_menu_help \"CSG Properties\" \"Defines appearence of current CSG geometry (color, visibility, transparency)\"\n",\
"set_menu_help \"STL Doctor\" \"Calls STL Doctor for preprocessing STL geometry files.\"\n",\
"set_menu_help \"STL Info\" \"Retrieves information about current STL geometry.\"\n",\
"\n",\
"set_menu_help \"Mesh\" \"Menu for mesh generation\"\n",\
"set_menu_help \"Generate Mesh\" \"Generates mesh from geometry, same as Button \\\"Generate Mesh\\\"\"\n",\
"set_menu_help \"Stop Meshing\" \"Terminates meshgeneration. It may take a while until meshing terminates, please be patient.\"\n",\
"set_menu_help \"Meshing Options\" \"Set options for mesh generation.\"\n",\
"set_menu_help \"Delete Mesh\" \"Deletes mesh. Not necessary before generation of new mesh.\"\n",\
"set_menu_help \"Delete Vol Mesh\" \"Deletes only volume mesh.\"\n",\
"set_menu_help \"Mesh Quality\" \"Computs element shape measures. Triangle angles are inner angles of all triangles (faces of tetrahedra). Tet angles are angles between faces of tetrahedra.\"\n",\
"set_menu_help \"Check Surface Mesh\" \"Checks consistency and overlap of surface mesh. Marks overlapping elements as bad elements, please enable visualization of bad elements in View->Mesh.\"\n",\
"set_menu_help \"Check Volume Mesh\" \"Checks conformity of volume mesh.\"\n",\
"set_menu_help \"Edit Boundary Conditions\" \"Open dialog for setting boundary condition numbers for individual faces.\"\n",\
"set_menu_help \"Analyze Geometry\" \"Perform only first step in mesh generation. Action depends on geometry type, e.g. generates charts for STL mesh, find vertices in CSG geometries.\"\n",\
"set_menu_help \"Mesh Edges\" \"Meshes edges\"\n",\
"set_menu_help \"Mesh Surface\" \"Generates surface mesh. Includes already surface optimization for some geomtry types.\"\n",\
"set_menu_help \"Optimize Surface\" \"Optimizes surface mesh.\"\n",\
"set_menu_help \"Surface Optim. Step\" \"Performs a specific surface optimiztion step. Mesh smoothing moves nodes. edge swapping swaps the diagonal of a quadrilateral built by two triangles, criterion either by number of nodes, or anlges. Combine points eliminates triangles by combining points (in the center of gravity).\"\n",\
"set_menu_help \"Mesh Volume\" \"Performs volume meshing. Algorithm is a combination of Delaunay and Rule-based Advancing Front\"\n",\
"set_menu_help \"Optimize Volume\" \"Performs additional volume optimization steps\"\n",\
"set_menu_help \"Smooth Opt Volume\" \"Performs optimization steps by smoothing iterations\"\n",\
"set_menu_help \"Smooth Opt Volume Jacobian\" \"Volume optimization by smoothing iterations. Criterion is optimization of Jacobi determinants. This optimization step is also available for 10-node tetrahedra.\"\n",\
"\n",\
"set_menu_help \"View\" \"Sets viewing options\"\n",\
"set_menu_help \"Zoom all\" \"Zooms scene to show whole object\"\n",\
"set_menu_help \"Center\" \"Defines center of rotation\"\n",\
"set_menu_help \"Viewing Options\" \"Sets viewing options for geometry, mesh, lighting\"\n",\
"set_menu_help \"Clipping Plane\" \"Introduces clipping plane. The clipping plane is defined by the normal vector, and a scaled offset. Clipping of performed by OpenGl rendering\"\n",\
"set_menu_help \"Quality Plot\" \"Shows the element quality distribution histogram. Measure is volume scaled by edge-length to the third. Optimal elements have measure 1.\"\n",\
"set_menu_help \"Sensitve Help\" \"Shows this help window\"\n",\
"\n",\
"set_menu_help \"Mesh-size\" \"Manipulations of existing mesh\"\n",\
"set_menu_help \"Refine uniform\" \"Refines mesh by splitting elements into eight childs (algorithm of J. Bey)\"\n",\
"set_menu_help \"Second Order\" \"Converts 4 node elements to 10 node elements. Edge-midpoitns are projected to the geometry.\"\n",\
"set_menu_help \"Refinement Dialog\" \"Controls local mesh refinement\"\n",\
"set_menu_help \"Load Meshsize\" \"Loads mesh-size file for local mesh refinement.\"\n",\
"set_menu_help \"MS from Surf Mesh\" \"Defines mesh-size by the surface mesh.\"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set f .options_dlg.nb.nbframe.general\n",\
"set_control_help $f.fine \"Controls relative mesh size.\\nThis control affects other mesh-size controls in common\"\n",\
"set_control_help $f.first \"First step in mesh generation. Usually, meshing  starts from \\\"analyze geometry\\\". If the surface mesh is already available \\\"First step\\\" should be set to \\\"mesh volume\\\"\"\n",\
"set_control_help $f.last \"Last step in mesh generation. If only the surface mesh is required, please set \\\"Last Step\\\" to \\\"Optimize Surface\\\"\"\n",\
"\n",\
"set_control_help .bubar.surfm \"Start mesh generation\"\n",\
"set_control_help .bubar.stopm \"Start mesh generation\"\n",\
"\n",\
"proc help_item { helptext } {p\n",\
"    puts $helptext\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc show_help { } {\n",\
"    \n",\
"    set w .help\n",\
"    \n",\
"    if {[winfo exists .help] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconif $w\n",\
"	focus $w \n",\
"    } {\n",\
"	\n",\
"	toplevel $w\n",\
"	\n",\
"	frame $w.buttons\n",\
"	pack $w.buttons -side bottom -fill x -pady 2m\n",\
"	button $w.buttons.done -text Done -command \"destroy $w\"\n",\
"	pack $w.buttons.done  -side left -expand 1\n",\
"\n",\
"	text $w.text -yscrollcommand \"$w.scroll set\" -setgrid true \\\n",\
"	    -width 60 -height 24 -wrap word\n",\
"	scrollbar $w.scroll -command \"$w.text yview\"\n",\
"	pack $w.scroll -side right -fill y\n",\
"	pack $w.text -expand yes -fill both\n",\
"\n",\
"    }\n",\
"    $w.text configure -state normal\n",\
"    $w.text delete 1.0 end\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set bold \"-background #43ce80 -relief raised -borderwidth 1\"\n",\
"set normal \"-background {} -relief flat\"\n",\
"\n",\
"\n",\
"proc help_main { } {\n",\
"\n",\
"    show_help;\n",\
"    set w .help\n",\
"    global bold\n",\
"    global normal\n",\
"\n",\
"\n",\
"    \n",\
"    $w.text insert 0.0 \\\n",\
"	{NETGEN Help}\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{1. General} d1\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{2. Menu items } d2\n",\
"    $w.text insert end \\n\\n\n",\
"\n",\
"    foreach tag {d1 d2} {\n",\
"	$w.text tag bind $tag <Any-Enter> \"$w.text tag configure $tag $bold\"\n",\
"	$w.text tag bind $tag <Any-Leave> \"$w.text tag configure $tag $normal\"\n",\
"    }\n",\
"    \n",\
"    $w.text tag bind d1 <1> { puts \"general\"; help_general }\n",\
"    $w.text tag bind d2 <1> { help_menus }\n",\
"\n",\
"    $w.text configure -state disabled\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc help_general { } {\n",\
"\n",\
"    show_help;\n",\
"    set w .help\n",\
"    global bold\n",\
"    global normal\n",\
"\n",\
"    puts \"general called\"\n",\
"\n",\
"    $w.text insert 0.0 \\\n",\
"	{NETGEN is an automatic three dimensional tetrahedral mesh generation system. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STEP or STL file format. NETGEN contains modules for mesh optimization and hierarchical mesh refinement.}\n",\
"\n",\
"    $w.text configure -state disabled\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc help_menus { } {\n",\
"\n",\
"    show_help;\n",\
"    set w .help\n",\
"    global bold\n",\
"    global normal\n",\
"\n",\
"\n",\
"    $w.text insert 0.0 \\\n",\
"	{The NETGEN Menu items are}\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{1. File} d1\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{2. Geometry } d2\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{3. Mesh } d3\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{4. View } d4\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{5. Mesh-size } d5\n",\
"    $w.text insert end \\n\\n\n",\
"    $w.text insert end \\\n",\
"	{6. STL } d6\n",\
"\n",\
"    foreach tag {d1 d2 d3 d4 d5 d6} {\n",\
"	$w.text tag bind $tag <Any-Enter> \"$w.text tag configure $tag $bold\"\n",\
"	$w.text tag bind $tag <Any-Leave> \"$w.text tag configure $tag $normal\"\n",\
"    }\n",\
"    \n",\
"    $w.text tag bind d1 <1> {puts \"File menu\"}\n",\
"    $w.text tag bind d2 <1> {puts \"Geometry menu\"}\n",\
"    $w.text tag bind d3 <1> {puts \"Mesh menu\"}\n",\
"    $w.text tag bind d4 <1> {puts \"View menu\"}\n",\
"    $w.text tag bind d5 <1> {puts \"Mesh-size menu\"}\n",\
"    $w.text tag bind d6 <1> {puts \"STL menu\"}\n",\
"\n",\
"    $w.text configure -state disabled\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"Ng_Vis_Set parameters\n",\
"\n",\
"\n",\
"\n",\
"set viscnt 0\n",\
"proc snapshottimer { } {\n",\
"    after 2000 { snapshottimer }\n",\
"\n",\
"    global viscnt\n",\
"    set viscnt [expr $viscnt+1]\n",\
"    set s1 0000$viscnt\n",\
"    set cnt [string range $s1 [expr [string length $s1]-4] end]\n",\
"    set filename \"p$cnt.jpg\"\n",\
"}\n",\
"snapshottimer\n",\
"\n",\
"\n",\
"\n",\
"proc redrawtimer { } {\n",\
"    global visoptions.autoredraw\n",\
"    global visoptions.autoredrawtime\n",\
"\n",\
"    set delay [expr int(${visoptions.autoredrawtime}*1000)]\n",\
"    if { ${visoptions.autoredraw} == 1 } { redraw; }\n",\
"    after $delay { redrawtimer } \n",\
"}\n",\
"redrawtimer\n",\
"\n",\
"\n",\
"set perstarttime [clock clicks -millisecond]\n",\
"proc redrawperiodic { } {\n",\
"    global visoptions.redrawperiodic\n",\
"    global perstarttime\n",\
"    set curtime [clock clicks -millisecond]\n",\
"    Ng_Vis_Set time [expr ($curtime - $perstarttime) / 5]\n",\
"    redraw\n",\
"    if { ${visoptions.redrawperiodic} == 1 } { after 30 { redrawperiodic } };\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc addplotline { identifier datax datay plotinfo {color black}} {\n",\
"    set c $identifier.c\n",\
"\n",\
"    set xstart [lindex $plotinfo 0]\n",\
"    set ystart [lindex $plotinfo 1]\n",\
"    set xmin [lindex $plotinfo 2]\n",\
"    set ymin [lindex $plotinfo 3]\n",\
"    set unitx [lindex $plotinfo 4]\n",\
"    set unity [lindex $plotinfo 5]\n",\
"    \n",\
"    \n",\
"    \n",\
"    set latestx [expr ([lindex $datax 0]-$xmin)*$unitx + $xstart]\n",\
"    set latesty [expr ([lindex $datay 0]-$ymin)*$unity + $ystart]\n",\
"    \n",\
"    for {set i 1} {$i < [llength $datax]} {incr i} {\n",\
"	set xpos [expr ([lindex $datax $i]-$xmin)*$unitx + $xstart]\n",\
"	set ypos [expr ([lindex $datay $i]-$ymin)*$unity + $ystart]\n",\
"	$c create line $latestx $latesty $xpos $ypos -width 1 -fill $color\n",\
"	set latestx $xpos\n",\
"	set latesty $ypos\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc createlineplot { width height identifier title xmin xmax ymin ymax plotinfo} {\n",\
"    set thiswidth $width\n",\
"    set thisheight $height\n",\
"    if { $thiswidth < 275 } { set thiswidth 275 }\n",\
"    if { $thisheight < 225 } { seth thisheight 225 }\n",\
"\n",\
"    set w $identifier\n",\
"\n",\
"    if {[winfo exists $w] == 1} {\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"	\n",\
"	set c $w.c\n",\
"	\n",\
"	canvas $c -relief raised -width $thiswidth -height $thisheight\n",\
"	pack $w.c -side top -fill x\n",\
"	\n",\
"	set titleFont {Helvetica 18}\n",\
"	set smallFont {Helvetica 12}\n",\
"\n",\
"	set xstart 100\n",\
"	set xend [expr $thiswidth-75]\n",\
"	set ystart [expr $thisheight-75]\n",\
"	set yend 75\n",\
"\n",\
"	$c create line $xstart $ystart $xstart $yend -width 2\n",\
"	$c create line $xstart $ystart $xend $ystart -width 2\n",\
"\n",\
"	\n",\
"\n",\
"	\n",\
"	set unitx [expr double($xend-$xstart)/($xmax-$xmin)]\n",\
"	set unity [expr double($yend-$ystart)/($ymax-$ymin)]\n",\
"\n",\
"	for {set i 0} {$i <= 1} {set i [expr $i+0.2]} {\n",\
"	    $c create line [expr $xstart+$i*($xend-$xstart)] [expr $ystart] [expr $xstart+$i*($xend-$xstart)] [expr $ystart+5] -width 2\n",\
"	    $c create text [expr $xstart+$i*($xend-$xstart)] [expr $ystart+7] -anchor n -font $smallFont \\\n",\
"		-text [format \"%.3g\" [expr $xmin+$i*($xmax-$xmin)]]\n",\
"	    $c create line [expr $xstart] [expr $ystart+$i*($yend-$ystart)] [expr $xstart-7] [expr $ystart+$i*($yend-$ystart)] -width 2\n",\
"	    $c create text [expr $xstart-9] [expr $ystart+$i*($yend-$ystart)] -anchor e -font $smallFont \\\n",\
"		-text [format \"%.3g\" [expr $ymin+$i*($ymax-$ymin)]]\n",\
"	}\n",\
"\n",\
"	upvar $plotinfo ploti\n",\
"\n",\
"	set ploti \"$xstart $ystart $xmin $ymin $unitx $unity\"\n",\
"\n",\
"		\n",\
"	button $w.close -text \"Close\" -command \"destroy $w\"\n",\
"	pack $w.close\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w $title\n",\
"	focus $w\n",\
"\n",\
"    }\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"proc getlineplotdata { datax datay xmini xmaxi ymini ymaxi} {\n",\
"\n",\
"    upvar $datax datx\n",\
"    upvar $datay daty\n",\
"    \n",\
"    upvar $xmini xmin\n",\
"    upvar $xmaxi xmax\n",\
"    upvar $ymini ymin\n",\
"    upvar $ymaxi ymax\n",\
"\n",\
"    global visoptions.lineplotusingx\n",\
"    global visoptions.lineplotusingy\n",\
"    global visoptions.lineplotsource\n",\
"    global visoptions.lineplotfile\n",\
"\n",\
"    set datx \"\"\n",\
"    set daty \"\"\n",\
"\n",\
"    set xmin 1e20\n",\
"    set xmax -1e20\n",\
"    set ymin 1e20\n",\
"    set ymax -1e20\n",\
"    \n",\
"\n",\
"    if {${visoptions.lineplotsource} == \"file\"} {\n",\
"	set fileId [open ${visoptions.lineplotfile} r]\n",\
"	set line \"\"\n",\
"	\n",\
"	while {[gets $fileId line] >= 0} {\n",\
"	    if { [string index [lindex $line 0] 0] != \"\\#\" } {\n",\
"		if { ${visoptions.lineplotusingx} < [llength $line] } {\n",\
"		    lappend datx [lindex $line ${visoptions.lineplotusingx}]\n",\
"		    \n",\
"		    if { [lindex $datx end] < $xmin } {set xmin [lindex $datx end]}\n",\
"		    if { [lindex $datx end] > $xmax } {set xmax [lindex $datx end]}\n",\
"		} {\n",\
"		    lappend datx 0\n",\
"		}\n",\
"		if { ${visoptions.lineplotusingy} < [llength $line] } {\n",\
"		    lappend daty [lindex $line ${visoptions.lineplotusingy}]\n",\
"\n",\
"		    if { [lindex $daty end] < $ymin } {set ymin [lindex $daty end]}\n",\
"		    if { [lindex $daty end] > $ymax } {set ymax [lindex $daty end]}\n",\
"		} {\n",\
"		    lappend daty 0\n",\
"		}\n",\
"	    }\n",\
"	    \n",\
"	}\n",\
"	close $fileId\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"proc lineplotdialog { } {\n",\
"\n",\
"    set w .lineplot_dlg\n",\
"    \n",\
"    if {[winfo exists .lineplot_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	\n",\
"	toplevel $w\n",\
"	\n",\
"	frame $w.filesettings -relief  groove -borderwidth 3\n",\
"	frame $w.filesettings.title\n",\
"	radiobutton $w.filesettings.title.choose -variable visoptions.lineplotsource \\\n",\
"	    -value file -text \"Data from File\"\n",\
"\n",\
"	pack $w.filesettings.title.choose -side left\n",\
"\n",\
"	pack $w.filesettings.title\n",\
"\n",\
"	\n",\
"	global visoptions.lineplotselectedeval\n",\
"	global visoptions.lineplotfile\n",\
"	global visoptions.evaluatefilenames\n",\
"	global visoptions.evaluatefiledescriptions\n",\
"\n",\
"	set evdata [NGS_GetData evaluatefiles]\n",\
"	set visoptions.evaluatefilenames none\n",\
"	set visoptions.evaluatefiledescriptions none\n",\
"	for {set i 0} {[expr $i+1] < [llength $evdata]} {incr i 2} {\n",\
"	    lappend visoptions.evaluatefilenames [lindex $evdata $i]\n",\
"	    lappend visoptions.evaluatefiledescriptions [lindex $evdata [expr $i+1]]	    \n",\
"	}\n",\
"	\n",\
"\n",\
"	tixOptionMenu $w.filesettings.latestevals -label \"Use Evaluate Results: \" \\\n",\
"	    -options {\n",\
"		label.width  25\n",\
"		label.anchor e\n",\
"		menubutton.width 40\n",\
"	    } \n",\
"	\n",\
"	for {set i 0} {$i < [llength ${visoptions.evaluatefilenames}]} {incr i} {\n",\
"	    $w.filesettings.latestevals add command $i \\\n",\
"		-label \"[lindex ${visoptions.evaluatefiledescriptions} $i] ([lindex ${visoptions.evaluatefilenames} $i])\"\n",\
"	}\n",\
"	$w.filesettings.latestevals config -variable visoptions.lineplotselectedeval\n",\
"\n",\
"	pack $w.filesettings.latestevals\n",\
"	\n",\
"	frame $w.filesettings.sfn\n",\
"\n",\
"	button $w.filesettings.sfn.bb -text \"Browse\" \\\n",\
"	    -command { set visoptions.lineplotfile [tk_getOpenFile] }\n",\
"\n",\
"	\n",\
"	entry $w.filesettings.sfn.fn -width 50 -relief sunken \\\n",\
"	    -textvariable visoptions.lineplotfile\n",\
"\n",\
"	pack $w.filesettings.sfn.bb $w.filesettings.sfn.fn -side left\n",\
"\n",\
"	pack $w.filesettings.sfn\n",\
"\n",\
"	button $w.filesettings.refresh -text \"Refresh\" -command {\n",\
"	    if { ${visoptions.lineplotselectedeval} != 0} {\n",\
"		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]\n",\
"	    }\n",\
"	    \n",\
"	    set saveusingx ${visoptions.lineplotusingx}\n",\
"	    set saveusingy ${visoptions.lineplotusingy}\n",\
"	    \n",\
"\n",\
"	    \n",\
"	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"		${visoptions.lineplotxcoordselector} delete $i\n",\
"	    }\n",\
"	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"		${visoptions.lineplotycoordselector} delete $i\n",\
"	    }\n",\
"\n",\
"	    \n",\
"	    set fileId [open ${visoptions.lineplotfile} r]\n",\
"	    set line \"\"\n",\
"	    gets $fileId line\n",\
"	    close $fileId\n",\
"	    if { [lindex $line 0] == \"\\#nglineplotinfo\" } {\n",\
"		set visoptions.lineplotdatadescr [lrange $line 1 end]	\n",\
"	    } {\n",\
"		set visoptions.lineplotdatadescr \"\"\n",\
"		for { set i 0 } { $i < [llength $line] } { incr i } {\n",\
"		    lappend visoptions.lineplotdatadescr \"data[expr $i+1]\"\n",\
"		}\n",\
"	    }\n",\
"\n",\
"	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"		${visoptions.lineplotxcoordselector} add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]\n",\
"	    }\n",\
"	    for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"		${visoptions.lineplotycoordselector} add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]\n",\
"	    }\n",\
"\n",\
"	    if { $saveusingx < [llength ${visoptions.lineplotdatadescr}] } {\n",\
"		set visoptions.lineplotusingx $saveusingx\n",\
"	    } {\n",\
"		set visoptions.lineplotusingx 0\n",\
"	    }\n",\
"	    if { $saveusingy < [llength ${visoptions.lineplotdatadescr}] } {\n",\
"		set visoptions.lineplotusingy $saveusingy\n",\
"	    } {\n",\
"		set visoptions.lineplotusingy 1\n",\
"	    }	    \n",\
"	}\n",\
"\n",\
"	pack $w.filesettings.refresh\n",\
"\n",\
"\n",\
"	frame $w.filesettings.using\n",\
"\n",\
"	global visoptions.lineplotdatadescr\n",\
"	\n",\
"	tixOptionMenu $w.filesettings.using.xco -label \"X-Coord:\"\\\n",\
"	    -options {\n",\
"		label.width  8\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"	for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"	    $w.filesettings.using.xco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]\n",\
"	}\n",\
"	$w.filesettings.using.xco config -variable visoptions.lineplotusingx\n",\
"	\n",\
"	tixOptionMenu $w.filesettings.using.yco -label \"Y-Coord:\"\\\n",\
"	    -options {\n",\
"		label.width  8\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"	for { set i 0 } { $i < [llength ${visoptions.lineplotdatadescr}] } { incr i } {\n",\
"	    $w.filesettings.using.yco add command $i -label [lindex ${visoptions.lineplotdatadescr} $i]\n",\
"	}\n",\
"	$w.filesettings.using.yco config -variable visoptions.lineplotusingy\n",\
"\n",\
"	global visoptions.lineplotxcoordselector\n",\
"	global visoptions.lineplotycoordselector\n",\
"	set visoptions.lineplotxcoordselector $w.filesettings.using.xco\n",\
"	set visoptions.lineplotycoordselector $w.filesettings.using.yco\n",\
"	\n",\
"\n",\
"	pack $w.filesettings.using.xco $w.filesettings.using.yco -side left\n",\
"	pack $w.filesettings.using\n",\
"\n",\
"	pack $w.filesettings -fill x -ipady 3\n",\
"	\n",\
"	frame $w.settings -relief  groove -borderwidth 3\n",\
"	label $w.settings.title -text \"\\nSettings\\n\"\n",\
"	pack $w.settings.title\n",\
"\n",\
"	frame $w.settings.minmax \n",\
"	checkbutton $w.settings.minmax.autoscale -text \"Autoscale\" -variable visoptions.lineplotautoscale\n",\
"	tixControl $w.settings.minmax.xmin -label \"Min. x: \" \\\n",\
"	    -integer false -variable visoptions.lineplotxmin \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 8\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	tixControl $w.settings.minmax.xmax -label \"Max. x: \" \\\n",\
"	    -integer false -variable visoptions.lineplotxmax \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 8\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	tixControl $w.settings.minmax.ymin -label \"Min. y: \" \\\n",\
"	    -integer false -variable visoptions.lineplotymin \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 8\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	tixControl $w.settings.minmax.ymax -label \"Max. y: \" \\\n",\
"	    -integer false -variable visoptions.lineplotymax \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 8\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	    \n",\
"\n",\
"	pack $w.settings.minmax.autoscale $w.settings.minmax.xmin $w.settings.minmax.xmax \\\n",\
"	    $w.settings.minmax.ymin $w.settings.minmax.ymax -side left\n",\
"\n",\
"	pack $w.settings.minmax\n",\
"\n",\
"	\n",\
"	label $w.settings.empty1 -text \"\"\n",\
"	pack $w.settings.empty1\n",\
"\n",\
"	frame $w.settings.plotsize\n",\
"\n",\
"	tixControl $w.settings.plotsize.xsize -label \"Plotsize  x: \"\\\n",\
"	    -integer true -variable visoptions.lineplotsizex \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 13\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	tixControl $w.settings.plotsize.ysize -label \"y: \"\\\n",\
"	    -integer true -variable visoptions.lineplotsizey \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 3\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	pack $w.settings.plotsize.xsize $w.settings.plotsize.ysize -side left\n",\
"\n",\
"	pack $w.settings.plotsize\n",\
"\n",\
"	label $w.settings.empty2 -text \"\"\n",\
"	pack $w.settings.empty2\n",\
"	\n",\
"\n",\
"	tixOptionMenu $w.settings.color -label \"Linecolor: \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"	foreach step { red black blue green yellow } {\n",\
"	    $w.settings.color add command $step -label $step\n",\
"	}\n",\
"	$w.settings.color config -variable visoptions.lineplotcolor\n",\
"\n",\
"	pack $w.settings.color\n",\
"\n",\
"\n",\
"	pack $w.settings\n",\
"\n",\
"	set datax \"\"\n",\
"	set datay \"\"\n",\
"	set xmin 0\n",\
"	set xmax 0\n",\
"	set ymin 0\n",\
"	set ymax 0\n",\
"\n",\
"	frame $w.plots -relief  groove -borderwidth 3\n",\
"\n",\
"	tixOptionMenu $w.plots.selplot -label \"Selected Plot: \" \\\n",\
"	    -options {\n",\
"		label.width  19\n",\
"		label.anchor e\n",\
"		menubutton.width 15\n",\
"	    } \n",\
"	$w.plots.selplot add command none -label \"None\"\n",\
"\n",\
"	$w.plots.selplot config -variable visoptions.lineplotselected\n",\
"\n",\
"	global visoptions.lineplotselector\n",\
"	set visoptions.lineplotselector $w.plots.selplot\n",\
"\n",\
"	\n",\
"\n",\
"	button $w.plots.new -text \"Generate New Plot\" -command {\n",\
"	    if { ${visoptions.lineplotselectedeval} != 0} {\n",\
"		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]\n",\
"	    }\n",\
"\n",\
"	    getlineplotdata datax datay xmin xmax ymin ymax\n",\
"\n",\
"	    puts stdout \"xmin $xmin xmax $xmax ymin $ymin ymax $ymax\"\n",\
"\n",\
"	    global visoptions.lineplotautoscale\n",\
"\n",\
"	    if {! ${visoptions.lineplotautoscale}} {\n",\
"		puts \"using set min/max values\"\n",\
"		set xmin ${visoptions.lineplotxmin}\n",\
"		set xmax ${visoptions.lineplotxmax}\n",\
"		set ymin ${visoptions.lineplotymin}\n",\
"		set ymax ${visoptions.lineplotymax}\n",\
"	    }\n",\
"	    \n",\
"	    incr visoptions.lineplotcurrentnum\n",\
"	    \n",\
"	    set ident .newplot${visoptions.lineplotcurrentnum}\n",\
"	    set plotinfo \"\"\n",\
"	    \n",\
"	    createlineplot ${visoptions.lineplotsizex} ${visoptions.lineplotsizey} \\\n",\
"		$ident \"Lineplot ${visoptions.lineplotcurrentnum}\" \\\n",\
"		$xmin $xmax $ymin $ymax plotinfo\n",\
"	    \n",\
"	    lappend visoptions.lineplotinfos $plotinfo\n",\
"\n",\
"	    \n",\
"	    ${visoptions.lineplotselector} add command ${visoptions.lineplotcurrentnum} -label \"Lineplot ${visoptions.lineplotcurrentnum}\"\n",\
"\n",\
"	    addplotline $ident $datax $datay $plotinfo ${visoptions.lineplotcolor}\n",\
"	}\n",\
"	\n",\
"	button $w.plots.addto -text \"Add to Selected Plot\" -command {\n",\
"	    if { ${visoptions.lineplotselectedeval} != 0} {\n",\
"		set visoptions.lineplotfile [lindex ${visoptions.evaluatefilenames} ${visoptions.lineplotselectedeval}]\n",\
"	    }\n",\
"\n",\
"	    if { ${visoptions.lineplotselected} != \"none\" } {\n",\
"\n",\
"		getlineplotdata datax datay xmin xmax ymin ymax\n",\
"\n",\
"		set ident .newplot${visoptions.lineplotselected}\n",\
"		set plotinfo [lindex ${visoptions.lineplotinfos} ${visoptions.lineplotselected}]\n",\
"\n",\
"		addplotline $ident $datax $datay $plotinfo ${visoptions.lineplotcolor}\n",\
"	    }\n",\
"	}\n",\
"	    \n",\
"	\n",\
"	\n",\
"	pack $w.plots.new $w.plots.addto $w.plots.selplot\n",\
"	\n",\
"	\n",\
"	pack $w.plots -fill x -ipady 3\n",\
"	\n",\
"\n",\
"\n",\
"	button $w.close -text \"Close\" -command \"destroy $w\"\n",\
"	pack $w.close\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +200+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"2D Lineplots\"\n",\
"\n",\
"	focus $w\n",\
"    \n",\
"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"set fieldlinesdialog_pop1 0\n",\
"\n",\
"\n",\
"proc fieldlinesdialog { } {\n",\
"    \n",\
"    set w .fieldlines_dlg\n",\
"\n",\
"    global fieldlinesdialog_pop1\n",\
"    set fieldlinesdialog_pop1 1\n",\
"\n",\
"    \n",\
"    if {[winfo exists .fieldlines_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	\n",\
"	toplevel $w\n",\
"\n",\
"	tixNoteBook $w.nb -ipadx 6 -ipady 6\n",\
"\n",\
"	$w.nb add draw -label \"Draw\"\n",\
"	$w.nb add settings -label \"Settings\"\n",\
"\n",\
"	pack $w.nb -expand yes -fill both -padx 5 -pady 5 -side top\n",\
"\n",\
"\n",\
"	\n",\
"	set f [$w.nb subwidget draw]\n",\
"	\n",\
"\n",\
"	frame $f.general\n",\
"\n",\
"	\n",\
"	checkbutton $f.general.enable -text \"Enable Fieldlines\" \\\n",\
"	    -variable visoptions.drawfieldlines \\\n",\
"	    -command { \n",\
"								Ng_Vis_Set parameters; \n",\
"		redraw \n",\
"	    }\n",\
"	tixControl $f.general.num -label \"Num: \" -integer true \\\n",\
"	    -variable visoptions.numfieldlines \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 12\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	\n",\
"	pack $f.general.enable $f.general.num -side left\n",\
"\n",\
"	pack $f.general\n",\
"\n",\
"	label $f.labe0 -text \" \"\n",\
"\n",\
"	pack $f.labe0\n",\
"\n",\
"	frame $f.general1\n",\
"	\n",\
"	checkbutton $f.general1.randomstart -text \"Field dependent density    \" \\\n",\
"	    -variable visoptions.fieldlinesrandomstart \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw}\n",\
"\n",\
"	\n",\
"	checkbutton $f.general1.redrawperiodic -text \"Animate periodic\" \\\n",\
"	    -variable visoptions.redrawperiodic \\\n",\
"	    -command { \n",\
"		redrawperiodic\n",\
"		Ng_Vis_Set parameters; \n",\
"		redraw \n",\
"	    }\n",\
"\n",\
"	pack $f.general1.randomstart $f.general1.redrawperiodic -side left\n",\
"\n",\
"	pack $f.general1\n",\
"\n",\
"\n",\
"\n",\
"	label $f.lab0 -text \" \"\n",\
"\n",\
"	pack $f.lab0\n",\
"\n",\
"\n",\
"	\n",\
"	tixOptionMenu $f.vecfun -label \"Vector Function: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"	$f.vecfun add command none -label None \n",\
"	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {\n",\
"	    set fname [Ng_Vis_Field getfieldname $i]\n",\
"	    set fcomp [Ng_Vis_Field getfieldcomponents $i]\n",\
"	    set iscomplex [Ng_Vis_Field iscomplex $i]\n",\
"	    set sdim [Ng_Vis_Field getdimension]\n",\
"	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }\n",\
"	    if { ($fcomp == $sdim) || ($fcomp == 3) } {\n",\
"		$f.vecfun add command $fname -label $fname\n",\
"	    } \n",\
"	}\n",\
"	$f.vecfun configure -variable visoptions.fieldlinesvecfunction\n",\
"	$f.vecfun configure -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	pack $f.vecfun\n",\
"	\n",\
"\n",\
"\n",\
"	label $f.lab00 -text \" \"\n",\
"\n",\
"	pack $f.lab00\n",\
"\n",\
"	frame $f.phasesettings\n",\
"\n",\
"	checkbutton $f.phasesettings.onephase -text \"Fix Phase\" -variable visoptions.fieldlinesonlyonephase\n",\
"	scale $f.phasesettings.phase -orient horizontal -length 300 -from 0 -to 360 \\\n",\
"	    -label \"phi\" \\\n",\
"	    -resolution 1 \\\n",\
"	    -variable visoptions.fieldlinesphase \\\n",\
"	    -command { popupcheckredraw3 fieldlinesdialog_pop1 }\n",\
"\n",\
"	pack $f.phasesettings.onephase $f.phasesettings.phase -side left\n",\
"	\n",\
"	pack $f.phasesettings\n",\
"\n",\
"\n",\
"\n",\
"	label $f.lab1 -text \" \"\n",\
"\n",\
"	pack $f.lab1\n",\
"\n",\
"\n",\
"	\n",\
"	frame $f.boxsettings -relief groove -borderwidth 3\n",\
"	frame $f.boxsettings.title\n",\
"	radiobutton $f.boxsettings.title.choose -variable visoptions.fieldlinesstartarea \\\n",\
"	    -value box -text \"Startpoints in Box\"\n",\
"\n",\
"	pack $f.boxsettings.title.choose -side left\n",\
"\n",\
"	pack $f.boxsettings.title\n",\
"\n",\
"	frame $f.boxsettings.points\n",\
"\n",\
"	label $f.boxsettings.points.lab2 -text \"Pmin\";\n",\
"	entry $f.boxsettings.points.ent1x -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap1x\n",\
"	entry $f.boxsettings.points.ent1y -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap1y\n",\
"	entry $f.boxsettings.points.ent1z -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap1z\n",\
"	label $f.boxsettings.points.lab3 -text \"   Pmax\";\n",\
"	entry $f.boxsettings.points.ent2x -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap2x\n",\
"	entry $f.boxsettings.points.ent2y -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap2y\n",\
"	entry $f.boxsettings.points.ent2z -width 8 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesstartareap2z\n",\
"	\n",\
"	pack $f.boxsettings.points\n",\
"	pack $f.boxsettings.points.lab2 $f.boxsettings.points.ent1x $f.boxsettings.points.ent1y $f.boxsettings.points.ent1z -side left\n",\
"	pack $f.boxsettings.points.lab3 $f.boxsettings.points.ent2x $f.boxsettings.points.ent2y $f.boxsettings.points.ent2z -side left\n",\
"\n",\
"	button $f.boxsettings.settobb -text \"Bounding Box\" -command {\n",\
"	    set bbox [Ng_MeshInfo bbox]\n",\
"	    set visoptions.fieldlinesstartareap1x [lindex $bbox 0]\n",\
"	    set visoptions.fieldlinesstartareap2x [lindex $bbox 1]\n",\
"	    set visoptions.fieldlinesstartareap1y [lindex $bbox 2]\n",\
"	    set visoptions.fieldlinesstartareap2y [lindex $bbox 3]\n",\
"	    set visoptions.fieldlinesstartareap1z [lindex $bbox 4]\n",\
"	    set visoptions.fieldlinesstartareap2z [lindex $bbox 5]\n",\
"	}\n",\
"\n",\
"	pack $f.boxsettings.settobb\n",\
"	    \n",\
"	\n",\
"	pack $f.boxsettings -fill x -ipady 3\n",\
"\n",\
"\n",\
"	frame $f.facesettings -relief groove -borderwidth 3\n",\
"	frame $f.facesettings.title\n",\
"	radiobutton $f.facesettings.title.choose -variable visoptions.fieldlinesstartarea \\\n",\
"	    -value face -text \"Startpoints on Face\"\n",\
"\n",\
"	pack $f.facesettings.title.choose -side left\n",\
"\n",\
"	pack $f.facesettings.title\n",\
"	\n",\
"	frame $f.facesettings.index\n",\
"	label $f.facesettings.index.lab -text \"face index:\"\n",\
"	label $f.facesettings.index.ent -text 1 -padx 4\n",\
"\n",\
"	pack $f.facesettings.index.lab $f.facesettings.index.ent -side left\n",\
"\n",\
"	pack $f.facesettings.index\n",\
"	\n",\
"	pack $f.facesettings -fill x -ipady 3\n",\
"\n",\
"\n",\
"	global visoptions.fieldlinesfilename\n",\
"\n",\
"	frame $f.filesettings -relief  groove -borderwidth 3\n",\
"	frame $f.filesettings.title\n",\
"	radiobutton $f.filesettings.title.choose -variable visoptions.fieldlinesstartarea \\\n",\
"	    -value file -text \"Startpoints from File\"\n",\
"\n",\
"	pack $f.filesettings.title.choose -side left\n",\
"\n",\
"	pack $f.filesettings.title\n",\
"\n",\
"	frame $f.filesettings.sfn\n",\
"\n",\
"	button $f.filesettings.sfn.bb -text \"Browse\" \\\n",\
"	    -command {\n",\
"		set types {\n",\
"		    { \"Netgen Fieldlines\" {.nef} }\n",\
"		}\n",\
"		set visoptions.fieldlinesfilename [tk_getOpenFile -filetypes $types -defaultextension \".nef\"]\n",\
"	    }\n",\
"\n",\
"	\n",\
"	entry $f.filesettings.sfn.fn -width 50 -relief sunken \\\n",\
"	    -textvariable visoptions.fieldlinesfilename\n",\
"\n",\
"	pack $f.filesettings.sfn.bb $f.filesettings.sfn.fn -side left\n",\
"\n",\
"	pack $f.filesettings.sfn\n",\
"\n",\
"	pack $f.filesettings -fill x -ipady 3\n",\
"	\n",\
"\n",\
"\n",\
"	\n",\
"	\n",\
"	set g [$w.nb subwidget settings]\n",\
"\n",\
"	frame $g.linesettings -relief groove -borderwidth 3\n",\
"	label $g.linesettings.title -text \"\\nLine Settings\\n\"\n",\
"	tixControl $g.linesettings.length -label \"rel. Length: \" -integer false \\\n",\
"	    -variable visoptions.fieldlineslength -min 0.00001 -max 10000 -step 0.1 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }\n",\
"\n",\
"	tixControl $g.linesettings.maxpoints -label \"max. Points: \" -integer true \\\n",\
"	    -variable visoptions.fieldlinesmaxpoints -min 0 -max 10000 -step 1 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }\n",\
"\n",\
"	tixControl $g.linesettings.thick -label \"rel. Thickness: \" -integer false \\\n",\
"	    -variable visoptions.fieldlinesthickness -min 1e-10 -max 0.5 -step 0.001 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }\n",\
"\n",\
"	pack $g.linesettings.title $g.linesettings.length $g.linesettings.maxpoints $g.linesettings.thick\n",\
"\n",\
"	pack $g.linesettings -fill x -ipady 3\n",\
"\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	global visoptions.fieldlinestolerance\n",\
"\n",\
"	frame $g.odesettings -relief groove -borderwidth 3\n",\
"	label $g.odesettings.title -text \"\\nODE Settings\\n\"\n",\
"	tixControl $g.odesettings.tol -label \"rel. Tolerance: \" -integer false \\\n",\
"	    -variable visoptions.fieldlinestolerance -min 0.00001 -max 1 -step 0.01 \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 25\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"\n",\
"	tixOptionMenu $g.odesettings.rktype -label \"RK-Type \" \\\n",\
"	    -options {\n",\
"		label.width  20\n",\
"		label.anchor e\n",\
"		menubutton.width 25\n",\
"	    }\n",\
"	$g.odesettings.rktype add command euler -label \"Euler, order 1\"\n",\
"	$g.odesettings.rktype add command eulercauchy -label \"Euler-Cauchy, order 2\"\n",\
"	$g.odesettings.rktype add command simpson -label \"Simpson, order 3\"\n",\
"	$g.odesettings.rktype add command crungekutta -label \"classical Runge-Kutta, order 4\"\n",\
"	$g.odesettings.rktype configure -variable visoptions.fieldlinesrktype\n",\
"	$g.odesettings.rktype configure -command { Ng_Vis_Set parameters; redraw }\n",\
"	\n",\
"	pack $g.odesettings.title $g.odesettings.tol $g.odesettings.rktype\n",\
"\n",\
"	pack $g.odesettings -fill x -ipady 3\n",\
"\n",\
"\n",\
"\n",\
"		\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu \n",\
"\n",\
"	button $w.bu.calc -text \"Build Fieldlines\" -command { \n",\
"	    if { ${visoptions.fieldlinesvecfunction} == \"none\" } {\n",\
"		bgerror \"Please select the vector function first!\"\n",\
"	    } {\n",\
"		set visoptions.drawfieldlines 1\n",\
"		Ng_Vis_Set parameters\n",\
"		Ng_BuildFieldLines\n",\
"		redraw \n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.bu.help -text \"Help\" -command {\n",\
"	    if {[winfo exists .fieldlines_help] == 1} {\n",\
"		wm withdraw .fieldlines_help\n",\
"		wm deiconify .fieldlines_help\n",\
"		focus .fieldlines_help\n",\
"	    } {\n",\
"		toplevel .fieldlines_help\n",\
"\n",\
"		tixScrolledText .fieldlines_help.ht -scrollbar y\n",\
"		set text [.fieldlines_help.ht subwidget text]\n",\
"\n",\
"		$text configure -setgrid true -wrap word \n",\
"\n",\
"		$text tag configure bold -font *-*-bold-*-*-*-*\n",\
"\n",\
"		\n",\
"		$text insert end \\\n",\
"		    \"Draw menu\\n \\n\" bold\n",\
"		$text insert end \\\n",\
"		    \"Enable Fieldlines\\n    To turn on and off the calculated fieldlines. (Has to be turned on to start the calculation)\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Num\\n    Number of fieldlines to calculate. (May not be used exactly.)\"\n",\
"		$text insert end \\\n",\
"		    \"Field dependent density\\n    There will be more fieldline startpoints where the field is stronger\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Animate periodic\\n    (for quasistationary fields) The fieldlines of the different phase angles are animated.\\n    ATTENTION: \\\"Fix Phase\\\" has to be turned off\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Vector Function\\n    The function fixing the direction of the lines\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Fix Phase\\n    (for quasistationary fields) Only calculate and draw fieldlines for one special phase angle.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Startpoints in Box\\n    Set the startpoints inside the box \\[Pmin1,Pmax1\\] x \\[Pmin2,Pmax2\\] x \\[Pmin3,Pmax3\\]\\n\"\n",\
"		$text insert end \\\n",\
"		    \"    With the button \\\"Bounding Box\\\" the whole bounding box of the geometry is selected.\\n\\n\" \n",\
"		$text insert end \\\n",\
"		    \"Startpoints on Face\\n    All startpoints will be set on one face. This face is selected by double-clicking with the mouse.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"Startpoints from File\\n    The startpoint information will be read from the selected file.\\n    The entries in the file can be as follows:\\n\"\n",\
"		$text insert end \\\n",\
"		    \"        point <x> <y> <z>\\n            set a (potential) startpoint\\n\"\n",\
"		$text insert end \\\n",\
"		    \"        line <x1> <y1> <z1> <x2> <y2> <z2> <n>\\n            set n (potential) startpoints on the line from (x1,y1,z1) to (x2,y2,z2)\\n\"\n",\
"		$text insert end \\\n",\
"		    \"        box <x1> <y1> <z1> <x2> <y2> <z2> <n>\\n            set n (potential) startpoints inside the box \\[x1,x2\\] x \\[y1,y2\\] x \\[z1,z2\\]\\n\"\n",\
"		$text insert end \\\n",\
"		    \"    ATTENTION: These are potential startpoints.\\n               The total number of startpoints will be bounded by the \\\"Num\\\"-parameter.\\n \\n \\n \\n\"\n",\
"		$text insert end \\\n",\
"		    \"Settings Menu\\n \\n\" bold\n",\
"		$text insert end \\\n",\
"		    \"rel. Length\\n    The maximal length of a fieldline relative to the diameter of the geometry.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"max. Points\\n    The maximum number of Runge-Kutta steps.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"rel. Thickness\\n    The thickness of the fieldlines relative to the diameter of the geometry.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"rel. Tolerance\\n    The tolerance for the step-length control of the Runge-Kutta method.\\n\\n\"\n",\
"		$text insert end \\\n",\
"		    \"RK-Type\\n    Which Runge-Kutta scheme to use\\n \\n \\n \\n\"\n",\
"		$text insert end \\\n",\
"		    \"Button \\\"Build Fieldlines\\\"\\n\" bold\n",\
"		$text insert end \\\n",\
"		    \"    Build the fieldlines.\"\n",\
"		\n",\
"\n",\
"		$text configure -state disabled\n",\
"\n",\
"		pack .fieldlines_help.ht -expand yes -fill both\n",\
"\n",\
"		wm withdraw .fieldlines_help\n",\
"		wm geom .fieldlines_help +300+200\n",\
"		wm deiconify .fieldlines_help\n",\
"		wm title .fieldlines_help \"Fieldlines Help\"\n",\
"		focus .fieldlines_help\n",\
"		\n",\
"	    }\n",\
"\n",\
"\n",\
"	}\n",\
"\n",\
"	button $w.bu.cancel -text \"Done\" -command \"destroy $w\"\n",\
"	pack $w.bu.calc $w.bu.help $w.bu.cancel -side left -expand yes\n",\
"	\n",\
"	\n",\
"	wm withdraw $w\n",\
"	wm geom $w +200+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Fieldlines\"\n",\
"		focus $w\n",\
"\n",\
"    }\n",\
"\n",\
"    global visoptions.fieldlinesstartface\n",\
"\n",\
"    \n",\
"    set f [$w.nb subwidget draw]\n",\
"    set visoptions.fieldlinesstartface [Ng_BCProp getactive]\n",\
"    $f.facesettings.index.ent configure -text ${visoptions.fieldlinesstartface}\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set visual_dialog_pop1 0\n",\
"set visual_dialog_pop2 0\n",\
"set visual_dialog_pop3 0\n",\
"set visual_dialog_pop4 0\n",\
"set visual_dialog_pop5 0\n",\
"set visual_dialog_pop6 0\n",\
"set visual_dialog_pop7 0\n",\
"\n",\
"proc visual_dialog { } {\n",\
"\n",\
"    set w .visoptions_dlg\n",\
"\n",\
"    \n",\
"    global visual_dialog_pop1\n",\
"    global visual_dialog_pop2\n",\
"    global visual_dialog_pop3\n",\
"    global visual_dialog_pop4\n",\
"    global visual_dialog_pop5\n",\
"    global visual_dialog_pop6\n",\
"    global visual_dialog_pop7\n",\
"    set visual_dialog_pop1 1\n",\
"    set visual_dialog_pop2 1\n",\
"    set visual_dialog_pop3 1\n",\
"    set visual_dialog_pop4 1\n",\
"    set visual_dialog_pop5 1\n",\
"    set visual_dialog_pop6 1\n",\
"    set visual_dialog_pop7 1\n",\
"   \n",\
"    if {[winfo exists .visoptions_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w\n",\
"    } {\n",\
"\n",\
"	toplevel $w\n",\
"\n",\
"	checkbutton $w.imaginary -text \"Imaginary Part\" \\\n",\
"	    -variable visoptions.imaginary \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	\n",\
"	frame $w.texframe\n",\
"\n",\
"	checkbutton $w.texframe.usetexture -text \"Use Textures (\" \\\n",\
"	    -variable visoptions.usetexture \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"	\n",\
"	checkbutton $w.texframe.lintexture -text \"Linear )\" \\\n",\
"	    -variable visoptions.lineartexture \\\n",\
"	    -command { Ng_Vis_Set parametersrange; redraw }\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	checkbutton $w.invcolor -text \"Inverse Color\" \\\n",\
"	    -variable visoptions.invcolor \\\n",\
"	    -command { Ng_Vis_Set parametersrange; redraw }\n",\
"	\n",\
"	checkbutton $w.redrawperiodic -text \"Animate periodic\" \\\n",\
"	    -variable visoptions.redrawperiodic \\\n",\
"	    -command { \n",\
"		redrawperiodic\n",\
"		Ng_Vis_Set parameters; \n",\
"		redraw \n",\
"	    }\n",\
"	\n",\
"\n",\
"	checkbutton $w.logscale -text \"Log Scale\" \\\n",\
"	    -variable visoptions.logscale \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"	\n",\
"	checkbutton $w.lineartexture -text \"Use Linear Texture\" \\\n",\
"	    -variable visoptions.lineartexture \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"	\n",\
"	scale $w.numcols -orient horizontal -length 100 -from 0 -to 50 \\\n",\
"	    -resolution 1   \\\n",\
"	    -variable  visoptions.numtexturecols \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop1 }\n",\
"\n",\
"	checkbutton $w.showclipsolution -text \"Draw Clipping Plane Solution\" \\\n",\
"	    -variable visoptions.showclipsolution \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	checkbutton $w.showsurfsolution -text \"Draw Surface Solution\" \\\n",\
"	    -variable visoptions.showsurfacesolution \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"	frame $w.grid -relief groove -borderwidth 3\n",\
"		scale $w.grid.size -orient horizontal -length 100 -from 1 -to 200 \\\n",\
"	    -label \"Grid\" \\\n",\
"	    -resolution 1    \\\n",\
"	    -variable  visoptions.gridsize \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop2 }\n",\
"    \n",\
"\n",\
"		scale $w.grid.xoffset -orient horizontal -length 80 -from 0 -to 1 \\\n",\
"	    -label \"x-Offset\" \\\n",\
"	    -resolution 0.05    \\\n",\
"	    -variable  visoptions.xoffset \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop3 }\n",\
"\n",\
"	scale $w.grid.yoffset -orient horizontal -length 80 -from 0 -to 1 \\\n",\
"	    -label \"y-Offset\" \\\n",\
"	    -resolution 0.05    \\\n",\
"	    -variable  visoptions.yoffset \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop4 }\n",\
"\n",\
"\n",\
"		pack $w.showsurfsolution\n",\
"	pack $w.grid -fill x -ipady 3\n",\
"	pack $w.grid.size $w.grid.xoffset $w.grid.yoffset -side left -expand yes\n",\
"\n",\
"\n",\
"\n",\
"	\n",\
"\n",\
"\n",\
"	frame $w.deform -relief groove -borderwidth 3\n",\
"	checkbutton $w.deform.cb -text \"Deformation\" \\\n",\
"	    -variable visoptions.deformation \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	tixControl $w.deform.sc1 -label \"Scale: \" -integer false \\\n",\
"	    -variable visoptions.scaledeform1 \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 7\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	scale $w.deform.sc2 -orient horizontal -length 100 -from 0 -to 1 \\\n",\
"	    -resolution 0.01    \\\n",\
"	    -variable  visoptions.scaledeform2 \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop5 }\n",\
"\n",\
"	pack $w.deform -fill x -ipady 2\n",\
"	pack $w.deform.cb $w.deform.sc1 $w.deform.sc2 -side left -expand yes\n",\
"	\n",\
"\n",\
"	frame $w.as -relief groove -borderwidth 3\n",\
"	checkbutton $w.as.autoscale -text \"Autoscale\" \\\n",\
"	    -variable visoptions.autoscale \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	tixControl $w.as.minval -label \"Min-value: \" -integer false \\\n",\
"	    -variable visoptions.mminval \\\n",\
"	    -command { Ng_Vis_Set parametersrange; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 12\n",\
"		label.anchor e\n",\
"	    }	\n",\
"	tixControl $w.as.maxval -label \"Max-value: \" -integer false \\\n",\
"	    -variable visoptions.mmaxval \\\n",\
"	    -command { Ng_Vis_Set parametersrange; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 12\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"	pack $w.as -fill x -ipady 3\n",\
"	pack $w.as.autoscale $w.as.minval $w.as.maxval -side left\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	frame $w.iso -relief groove -borderwidth 3\n",\
"	pack $w.iso -fill x -ipady 3\n",\
"\n",\
"	frame $w.iso.cb\n",\
"	pack $w.iso.cb -side left\n",\
"\n",\
"	checkbutton $w.iso.cb.isolines -text \"Iso-lines\" \\\n",\
"	    -variable visoptions.isolines \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"	pack $w.iso.cb.isolines -side top\n",\
"\n",\
"	checkbutton $w.iso.cb.isosurf -text \"Iso-Surface\" \\\n",\
"	    -variable visoptions.isosurf \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"	pack $w.iso.cb.isosurf -side top\n",\
"\n",\
"\n",\
"\n",\
"	scale $w.iso.numiso -orient horizontal -length 100 -from 2 -to 50 \\\n",\
"	    -label \"\" \\\n",\
"	    -resolution 1    \\\n",\
"	    -variable  visoptions.numiso \\\n",\
"	    -command { popupcheckredraw visual_dialog_pop6 }\n",\
"\n",\
"	pack $w.iso.numiso -side left\n",\
"\n",\
"\n",\
"\n",\
"	frame $w.iso.subdiv\n",\
"	radiobutton $w.iso.subdiv.zero -text \"0\" -variable visoptions.subdivisions -value 0 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	}\n",\
"	radiobutton $w.iso.subdiv.one -text \"1\" -variable visoptions.subdivisions -value 1 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	}\n",\
"	radiobutton $w.iso.subdiv.two -text \"2\" -variable visoptions.subdivisions -value 2 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	}\n",\
"	radiobutton $w.iso.subdiv.three -text \"3\" -variable visoptions.subdivisions -value 3 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	}\n",\
"	radiobutton $w.iso.subdiv.four -text \"4\" -variable visoptions.subdivisions -value 4 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	}\n",\
"	radiobutton $w.iso.subdiv.five -text \"5\" -variable visoptions.subdivisions -value 5 \\\n",\
"	    -command { \n",\
"				Ng_Vis_Set parameters; redraw;\n",\
"	    }\n",\
"\n",\
" 	label $w.iso.subdiv.text  -text \"subdivision\"\n",\
"\n",\
"	pack $w.iso.subdiv -side right -ipadx 10\n",\
"\n",\
"	pack $w.iso.subdiv.text -side top\n",\
"	pack $w.iso.subdiv.zero $w.iso.numiso -side left\n",\
"	pack $w.iso.subdiv.one $w.iso.numiso -side left\n",\
"	pack $w.iso.subdiv.two $w.iso.numiso -side left\n",\
"	pack $w.iso.subdiv.three $w.iso.numiso -side left\n",\
"	pack $w.iso.subdiv.four $w.iso.numiso -side left\n",\
"	pack $w.iso.subdiv.five $w.iso.numiso -side left\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	checkbutton $w.showcurves -text \"Show Curves\" \\\n",\
"	    -variable visoptions.drawpointcurves \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	pack $w.showcurves\n",\
"\n",\
"	frame $w.redraw -relief groove -borderwidth 3\n",\
"	checkbutton $w.redraw.auto -text \"Auto-redraw\" \\\n",\
"	    -variable visoptions.autoredraw \n",\
"\n",\
"	tixControl $w.redraw.val -label \" after (sec) \" -integer false \\\n",\
"	    -variable visoptions.autoredrawtime \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 0\n",\
"		label.anchor w\n",\
"	    }	\n",\
"\n",\
"	pack $w.redraw -fill x -ipady 3\n",\
"	pack $w.redraw.auto  $w.redraw.val -side left\n",\
"\n",\
"\n",\
"\n",\
"	tixControl $w.redraw.simtime -label \" Simulation Time (1e-6 s)\" -integer false \\\n",\
"	    -variable visoptions.simulationtime \\\n",\
"	    -command { \n",\
"		Ng_Vis_Set time ${visoptions.simulationtime}; \n",\
"		catch {NGS_Set time ${visoptions.simulationtime};}\n",\
"		redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 0\n",\
"		label.anchor w\n",\
"	    }	\n",\
"	pack $w.redraw.simtime -side left\n",\
"	\n",\
"\n",\
"\n",\
"\n",\
"	tixOptionMenu $w.clipsol -label \"Clipping Plane Sol: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"\n",\
"	set none 1\n",\
"	$w.clipsol add command none -label None\n",\
"	$w.clipsol add command scal -label \"Scalar Function\"\n",\
"	$w.clipsol add command vec -label \"Vector Function\"\n",\
"\n",\
"	$w.clipsol configure -variable visoptions.clipsolution\n",\
"	$w.clipsol configure -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	pack $w.clipsol\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	tixOptionMenu $w.scalfun -label \"Scalar Function: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"\n",\
"	tixOptionMenu $w.vecfun -label \"Vector Function: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"\n",\
"\n",\
"	$w.scalfun add command none -label None\n",\
"	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {\n",\
"	    set fname [Ng_Vis_Field getfieldname $i]\n",\
"	    set fcomp [Ng_Vis_Field getfieldcomponents $i]\n",\
"	    if { $fcomp == 1 } {\n",\
"		$w.scalfun add command $fname.1 -label $fname\n",\
"	    } {\n",\
"		for { set j 1 } { $j <= $fcomp } { incr j } {\n",\
"		    $w.scalfun add command $fname.$j -label \"$fname ($j)\"\n",\
"		}\n",\
"		$w.scalfun add command $fname.0 -label \"func ($fname)\"\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	$w.vecfun add command none -label None \n",\
"	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {\n",\
"	    set fname [Ng_Vis_Field getfieldname $i]\n",\
"	    set fcomp [Ng_Vis_Field getfieldcomponents $i]\n",\
"	    set iscomplex [Ng_Vis_Field iscomplex $i]\n",\
"	    set sdim [Ng_Vis_Field getdimension]\n",\
"	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }\n",\
"	    if { ($fcomp == $sdim) || ($fcomp == 3) } {\n",\
"		$w.vecfun add command $fname -label $fname\n",\
"	    } \n",\
"	}\n",\
"\n",\
"\n",\
"	$w.scalfun configure -variable visoptions.scalfunction \n",\
"	$w.scalfun configure -command { Ng_Vis_Set parameters; redraw }\n",\
"	$w.vecfun configure -variable visoptions.vecfunction\n",\
"	$w.vecfun configure -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	tixOptionMenu $w.evaluate -label \"Evaluate: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }	\n",\
"	$w.evaluate add command abs -label \"|.|\"\n",\
"	$w.evaluate add command abstens -label \"|tensor|\"\n",\
"	$w.evaluate add command mises -label \"Mises\"\n",\
"	$w.evaluate add command main  -label \"Main\"\n",\
"	$w.evaluate configure -variable visoptions.evaluate\n",\
"	$w.evaluate configure -command { \n",\
"	    Ng_Vis_Set parameters; \n",\
"	    redraw \n",\
"	}\n",\
"\n",\
"	pack $w.scalfun $w.vecfun $w.evaluate\n",\
"\n",\
"	tixControl $w.multidimcomp -label \"multidim-component: \" -integer true \\\n",\
"	    -variable visoptions.multidimcomponent -min 0 \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 18\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"\n",\
"	pack $w.multidimcomp\n",\
"\n",\
"	pack $w.imaginary $w.logscale $w.texframe $w.invcolor $w.redrawperiodic\n",\
"	pack $w.texframe.usetexture $w.texframe.lintexture -side left -expand yes\n",\
"	\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu  -pady 5\n",\
"\n",\
"	button $w.bu.showsol -text \"Show Solution\" -command { \n",\
"	    set selectvisual solution\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	}\n",\
"	button $w.bu.clipping -text \"Clipping\" -command { \n",\
"	    clippingdialog; \n",\
"	}\n",\
"	button $w.bu.fieldlines -text \"Fieldlines\" -command { \n",\
"	    fieldlinesdialog; \n",\
"	}\n",\
"\n",\
"	button $w.bu.lineplot -text \"2D Lineplot\" -command {\n",\
"	    lineplotdialog;\n",\
"	}\n",\
"\n",\
"	button $w.bu.done -text \"Close\" -command { \n",\
"	    destroy .visoptions_dlg\n",\
"	}\n",\
"\n",\
"	pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +100+100\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Visualization\"\n",\
"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc reset_visual_dialog { } {\n",\
"    \n",\
"    set w .visoptions_dlg\n",\
"    \n",\
"      if {[winfo exists .visoptions_dlg] == 1} {\n",\
"    \n",\
"    \n",\
"	  destroy $w.scalfun $w.vecfun $w.evaluate $w.multidimcomp\n",\
"	  destroy $w.imaginary $w.logscale $w.texframe.usetexture $w.texframe.lintexture\n",\
"          destroy $w.texframe\n",\
"          destroy $w.invcolor $w.redrawperiodic\n",\
"	  destroy $w.bu  -pady 5\n",\
"	  destroy $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes\n",\
"	  \n",\
"	  \n",\
"	  checkbutton $w.imaginary -text \"Imaginary Part\" \\\n",\
"	      -variable visoptions.imaginary \\\n",\
"	      -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	  frame $w.texframe\n",\
"\n",\
"	  checkbutton $w.texframe.usetexture -text \"Use Textures (\" \\\n",\
"	      -variable visoptions.usetexture \\\n",\
"	      -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"	  checkbutton $w.texframe.lintexture -text \"Linear )\" \\\n",\
"	      -variable visoptions.lineartexture \\\n",\
"	      -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	  \n",\
"	  \n",\
"	  checkbutton $w.invcolor -text \"Inverse Color\" \\\n",\
"	      -variable visoptions.invcolor \\\n",\
"	      -command { Ng_Vis_Set parameters; redraw }\n",\
"	  \n",\
"	  checkbutton $w.redrawperiodic -text \"Animate periodic\" \\\n",\
"	      -variable visoptions.redrawperiodic \\\n",\
"	      -command { \n",\
"		  redrawperiodic\n",\
"		  Ng_Vis_Set parameters; \n",\
"		  redraw \n",\
"	      }\n",\
"	  \n",\
"\n",\
"	checkbutton $w.logscale -text \"Log Scale\" \\\n",\
"	    -variable visoptions.logscale \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"\n",\
"	tixOptionMenu $w.scalfun -label \"Scalar Function: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"\n",\
"	tixOptionMenu $w.vecfun -label \"Vector Function: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }\n",\
"\n",\
"\n",\
"\n",\
"	$w.scalfun add command none -label None\n",\
"	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {\n",\
"	    set fname [Ng_Vis_Field getfieldname $i]\n",\
"	    set fcomp [Ng_Vis_Field getfieldcomponents $i]\n",\
"	    if { $fcomp == 1 } {\n",\
"		$w.scalfun add command $fname.1 -label $fname\n",\
"	    } {\n",\
"		for { set j 1 } { $j <= $fcomp } { incr j } {\n",\
"		    $w.scalfun add command $fname.$j -label \"$fname ($j)\"\n",\
"		}\n",\
"		$w.scalfun add command $fname.0 -label \"func ($fname)\"\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	$w.vecfun add command none -label None \n",\
"	for { set i 1 } { $i <= [Ng_Vis_Field getnfieldnames] } { incr i } {\n",\
"	    set fname [Ng_Vis_Field getfieldname $i]\n",\
"	    set fcomp [Ng_Vis_Field getfieldcomponents $i]\n",\
"	    set iscomplex [Ng_Vis_Field iscomplex $i]\n",\
"	    set sdim [Ng_Vis_Field getdimension]\n",\
"	    if { $iscomplex == 1 } { set fcomp [expr $fcomp / 2] }\n",\
"	    if { ($fcomp == $sdim) || ($fcomp == 3) } {\n",\
"		$w.vecfun add command $fname -label $fname\n",\
"	    } \n",\
"	}\n",\
"\n",\
"\n",\
"\n",\
"	$w.scalfun configure -variable visoptions.scalfunction \n",\
"	$w.scalfun configure -command { Ng_Vis_Set parameters; redraw }\n",\
"	$w.vecfun configure -variable visoptions.vecfunction\n",\
"	$w.vecfun configure -command { Ng_Vis_Set parameters; redraw }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	tixOptionMenu $w.evaluate -label \"Evaluate: \" \\\n",\
"	    -options {\n",\
"		label.width  18\n",\
"		label.anchor e\n",\
"		menubutton.width 12\n",\
"	    }	\n",\
"	$w.evaluate add command abs -label \"|.|\"\n",\
"	$w.evaluate add command abstens -label \"|tensor|\"\n",\
"	$w.evaluate add command mises -label \"Mises\"\n",\
"	$w.evaluate add command main  -label \"Main\"\n",\
"	$w.evaluate configure -variable visoptions.evaluate\n",\
"	$w.evaluate configure -command { \n",\
"	    Ng_Vis_Set parameters; \n",\
"	    redraw \n",\
"	}\n",\
"\n",\
"	pack $w.scalfun $w.vecfun $w.evaluate\n",\
"\n",\
"	tixControl $w.multidimcomp -label \"multidim-component: \" -integer true \\\n",\
"	    -variable visoptions.multidimcomponent -min 0 \\\n",\
"	    -command { Ng_Vis_Set parameters; redraw } \\\n",\
"	    -options {\n",\
"		entry.width 6\n",\
"		label.width 18\n",\
"		label.anchor e\n",\
"	    }	\n",\
"\n",\
"\n",\
"          pack $w.multidimcomp\n",\
"\n",\
"          pack $w.imaginary $w.logscale $w.texframe $w.invcolor $w.redrawperiodic\n",\
"          pack $w.texframe.usetexture $w.texframe.lintexture  -side left -expand yes\n",\
"\n",\
"\n",\
"	frame $w.bu\n",\
"	pack $w.bu  -pady 5\n",\
"\n",\
"	button $w.bu.showsol -text \"Show Solution\" -command { \n",\
"	    set selectvisual solution\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	}\n",\
"	button $w.bu.clipping -text \"Clipping\" -command { \n",\
"	    clippingdialog; \n",\
"	}\n",\
"	button $w.bu.fieldlines -text \"Fieldlines\" -command { \n",\
"	    fieldlinesdialog; \n",\
"	}\n",\
"\n",\
"	button $w.bu.lineplot -text \"2D Lineplot\" -command {\n",\
"	    lineplotdialog;\n",\
"	}\n",\
"\n",\
"	button $w.bu.done -text \"Close\" -command { \n",\
"	    destroy .visoptions_dlg\n",\
"	}\n",\
"\n",\
"	pack $w.bu.showsol $w.bu.clipping $w.bu.fieldlines $w.bu.lineplot $w.bu.done -side left -expand yes\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"\n",\
"\n",\
"      }\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"set sockets.serverport 0\n",\
"set sockets.serverhost \"localhost\"\n",\
"set sockets.serverlistbox 0\n",\
"set sockets.queuelistbox 0\n",\
"set sockets.currentjoblistbox 0\n",\
"set sockets.answerlistbox 0\n",\
"set sockets.myidlabel -1\n",\
"\n",\
"\n",\
"proc updateserverlist { } {\n",\
"    global sockets.serverlistbox\n",\
"    \n",\
"    set retval [Ng_Socket getserverlist]\n",\
"\n",\
"    ${sockets.serverlistbox} delete 0 end\n",\
"\n",\
"    for {set i 0} {$i < [llength $retval]} {incr i 3} {\n",\
"	${sockets.serverlistbox} insert end \\\n",\
"	    [format \"%-16s   %6i   %6i\" [lindex $retval $i] [lindex $retval [expr $i+1]] [lindex $retval [expr $i+2]]]\n",\
"    }\n",\
"}\n",\
"\n",\
"proc clientsocketdialog { } {\n",\
"    set w .clientsock_dlg\n",\
"    \n",\
"    if {[winfo exists .clientsock_dlg] == 1} {\n",\
"	wm withdraw $w\n",\
"	wm deiconify $w\n",\
"	focus $w \n",\
"    } {\n",\
"	toplevel $w\n",\
"\n",\
"	global sockets.serverhost\n",\
"	global sockets.serverport\n",\
"\n",\
"	frame $w.general\n",\
"	frame $w.host\n",\
"	label $w.host.lab -text \"Serverhost: \"\n",\
"	entry $w.host.name -width 30 -relief sunken -textvariable sockets.serverhost\n",\
"\n",\
"	pack $w.host.lab $w.host.name -side left\n",\
"	pack $w.host\n",\
"\n",\
"	frame $w.ports\n",\
"	label $w.ports.lab1 -text \"Serverport: \"\n",\
"	entry $w.ports.statport -width 6 -relief sunken -textvariable sockets.serverport\n",\
"	\n",\
"	pack $w.ports.lab1 $w.ports.statport -side left\n",\
"	pack $w.ports\n",\
"\n",\
"	frame $w.listboxes\n",\
"\n",\
"	frame $w.listboxes.choosesocketframe\n",\
"\n",\
"	tixScrolledListBox $w.listboxes.choosesocketframe.choosesocket -scrollbar auto\n",\
"\n",\
"	global sockets.serverlistbox\n",\
"\n",\
"	set sockets.serverlistbox [$w.listboxes.choosesocketframe.choosesocket subwidget listbox]\n",\
"\n",\
"	${sockets.serverlistbox} configure -width 35\n",\
"	${sockets.serverlistbox} configure -selectmode browse\n",\
"	${sockets.serverlistbox} configure -exportselection false\n",\
"\n",\
"	button $w.addserver -text \"Add ServerSocket\" -command {\n",\
"	    Ng_Socket addserver ${sockets.serverport} ${sockets.serverhost}\n",\
"	    updateserverlist\n",\
"	}\n",\
"	\n",\
"	pack $w.addserver\n",\
"\n",\
"	label $w.linefeed -text \"\\n\"\n",\
"	pack $w.linefeed\n",\
"	\n",\
"	frame $w.clientidframe\n",\
"	label $w.clientidframe.lab -text \"Client ID: \";\n",\
"	global sockets.myidlabel\n",\
"	entry $w.clientidframe.val -width 5 -relief sunken -textvariable sockets.myidlabel\n",\
"	button $w.clientidframe.but -text \"Set\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		Ng_Socket setid $opserver ${sockets.myidlabel}\n",\
"		updateserverlist\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	pack $w.clientidframe.lab $w.clientidframe.val $w.clientidframe.but -side left\n",\
"	pack $w.clientidframe\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	label $w.listboxes.choosesocketframe.chooselab -text [format \"\\n\\n%-16s    %6s  %6s                       \" Host Socket MyID ]\n",\
"	pack $w.listboxes.choosesocketframe.chooselab\n",\
"	pack $w.listboxes.choosesocketframe.choosesocket\n",\
"\n",\
"	frame $w.listboxes.choosesocketframe.serverbuttons\n",\
"\n",\
"	button $w.listboxes.choosesocketframe.serverbuttons.save -text \"Save\" -command {\n",\
"	    Ng_Socket saveserverlist\n",\
"	}\n",\
"\n",\
"	global sockets.serverlist\n",\
"	Ng_Socket loadserverlist\n",\
"	updateserverlist\n",\
"\n",\
"	button $w.listboxes.choosesocketframe.serverbuttons.delete -text \"Delete\" -command {\n",\
"	   set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		Ng_Socket deletesocket [lindex $opsel 0]\n",\
"		updateserverlist\n",\
"	    } \n",\
"	}\n",\
"	\n",\
"	pack $w.listboxes.choosesocketframe.serverbuttons.save $w.listboxes.choosesocketframe.serverbuttons.delete -side left\n",\
"	pack $w.listboxes.choosesocketframe.serverbuttons\n",\
"\n",\
"	frame $w.listboxes.statusframe\n",\
"\n",\
"	label $w.listboxes.statusframe.statuslabel1 -text \"\\n\\njobqueue\"\n",\
"\n",\
"	tixScrolledListBox $w.listboxes.statusframe.queuestatus -scrollbar auto\n",\
"\n",\
"	label $w.listboxes.statusframe.statuslabel2 -text \"\\ncurrent job\"\n",\
"\n",\
"	tixScrolledListBox $w.listboxes.statusframe.currentjobstatus -scrollbar auto\n",\
"\n",\
"	label $w.listboxes.statusframe.statuslabel3 -text \"\\nanswers\"\n",\
"\n",\
"	tixScrolledListBox $w.listboxes.statusframe.answers -scrollbar auto\n",\
"\n",\
"	global sockets.queuelistbox\n",\
"	global sockets.currentjoblistbox\n",\
"	global sockets.answerlistbox\n",\
"\n",\
"	set sockets.queuelistbox [$w.listboxes.statusframe.queuestatus subwidget listbox]\n",\
"	set sockets.currentjoblistbox [$w.listboxes.statusframe.currentjobstatus subwidget listbox]\n",\
"	set sockets.answerlistbox [$w.listboxes.statusframe.answers subwidget listbox]\n",\
"\n",\
"	${sockets.queuelistbox} configure -width 50\n",\
"	${sockets.queuelistbox} configure -height 5\n",\
"	${sockets.queuelistbox} configure -selectmode browse\n",\
"	${sockets.queuelistbox} configure -exportselection false\n",\
"	\n",\
"	${sockets.currentjoblistbox} configure -width 50\n",\
"	${sockets.currentjoblistbox} configure -height 1\n",\
"	${sockets.currentjoblistbox} configure -selectmode browse\n",\
"	${sockets.currentjoblistbox} configure -exportselection false\n",\
"\n",\
"	${sockets.answerlistbox} configure -width 50\n",\
"	${sockets.answerlistbox} configure -height 5\n",\
"	${sockets.answerlistbox} configure -selectmode browse\n",\
"	${sockets.answerlistbox} configure -exportselection false\n",\
"\n",\
"	button $w.listboxes.statusframe.updatebutton -text \"Update\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [Ng_Socket sendqueuestatus $opserver]\n",\
"\n",\
"		${sockets.queuelistbox} delete 0 end\n",\
"		\n",\
"		if {[lindex $retval 0] > 0} {\n",\
"		    ${sockets.queuelistbox} insert end [format \"Blocked for user %i\" [lindex $retval 0]]\n",\
"		} {\n",\
"		    ${sockets.queuelistbox} insert end \"Not blocked\"\n",\
"		}\n",\
"		\n",\
"		for {set i 2} {$i < [expr 2*[lindex $retval 1]+2]} {incr i 2} {\n",\
"		    ${sockets.queuelistbox} insert end [format \"client %i, command %s\" [lindex $retval $i] [lindex $retval [expr $i+1]]]\n",\
"		}\n",\
"		\n",\
"		${sockets.answerlistbox} delete 0 end\n",\
"		\n",\
"		for {set i [expr 2*[lindex $retval 1]+3]} {$i < [llength $retval]} {incr i 2} {\n",\
"		    ${sockets.answerlistbox} insert end [format \"client %i, command %s\" [lindex $retval $i] [lindex $retval [expr $i+1]]]\n",\
"		}\n",\
"\n",\
"		${sockets.currentjoblistbox} delete 0 end\n",\
"		set retval [Ng_Socket sendjobstatus $opserver]\n",\
"		if {[lindex $retval 0] != 0} {\n",\
"		    ${sockets.currentjoblistbox} insert end [format \"client %i, command %s: %s\" [lindex $retval 0] [lindex $retval 1] [lrange $retval 2 end]]\n",\
"		}\n",\
"		\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	pack $w.listboxes.statusframe.statuslabel1 $w.listboxes.statusframe.queuestatus \\\n",\
"	    $w.listboxes.statusframe.statuslabel2 $w.listboxes.statusframe.currentjobstatus \\\n",\
"	    $w.listboxes.statusframe.statuslabel3 $w.listboxes.statusframe.answers \\\n",\
"	    $w.listboxes.statusframe.updatebutton\n",\
"\n",\
"	pack $w.listboxes.choosesocketframe $w.listboxes.statusframe -side left\n",\
"	\n",\
"	pack $w.listboxes\n",\
"\n",\
"	label $w.lab1 -text \"\\n\"\n",\
"	pack $w.lab1\n",\
"\n",\
"\n",\
"	frame $w.buttons1\n",\
"	frame $w.buttons2\n",\
"	\n",\
"	button $w.buttons1.getid -text \"Get ID\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [Ng_Socket getid $opserver]\n",\
"		updateserverlist\n",\
"		set sockets.myidlabel $retval\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons1.killjob -text \"Kill Cur. Job\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		Ng_Socket killcurrentjob $opserver\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons2.sendmesh -text \"Send Mesh\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [Ng_Socket sendmesh $opserver]\n",\
"		set sockets.meshsent 1\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons2.sendpde -text \"Send PDE\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [NGS_Socket sendpdefile $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons2.solvepde -text \"Solve PDE\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [NGS_Socket solvepde $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons2.writesol -text \"Write Solution\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [NGS_Socket writesolution $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons2.sendsol -text \"Receive Solution\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [NGS_Socket sendsolution $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons1.blockserver -text \"Block Server\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [Ng_Socket blockserver $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	button $w.buttons1.unblockserver -text \"UnBlock Server\" -command {\n",\
"	    set opsel [${sockets.serverlistbox} curselection]\n",\
"	    if {[llength $opsel] > 0} {\n",\
"		set opserver [lindex $opsel 0]\n",\
"		set retval [Ng_Socket unblockserver $opserver]\n",\
"	    }\n",\
"	}\n",\
"\n",\
"	\n",\
"	pack $w.buttons1.getid $w.buttons1.blockserver $w.buttons1.unblockserver $w.buttons1.killjob -side left\n",\
"	pack $w.buttons2.sendmesh $w.buttons2.sendpde $w.buttons2.solvepde $w.buttons2.writesol $w.buttons2.sendsol -side left\n",\
"\n",\
"	pack $w.buttons1 $w.buttons2\n",\
"\n",\
"\n",\
"	wm withdraw $w\n",\
"	wm geom $w +200+200\n",\
"	wm deiconify $w\n",\
"	wm title $w \"Client Socket\"\n",\
"	focus .options_dlg\n",\
"\n",\
"    }\n",\
"    \n",\
"\n",\
"}\n",\
"\n",\
".ngmenu.special add command -label \"Client Socket\" \\\n",\
"    -command { clientsocketdialog }\n",\
"\n",\
"\n",\
"\n",\
"set entities [ ]\n",\
"\n",\
"\n",\
"proc acisdialogbuildtree {} {\n",\
"\n",\
"    global entities\n",\
"\n",\
"    set w .acis_dlg\n",\
"    set hlist [$w.mtre subwidget hlist]\n",\
"    \n",\
"    set entities [Ng_ACISCommand getentities]\n",\
"\n",\
"    set nrentities [expr [llength $entities]]\n",\
"\n",\
"    if {$nrentities != 0} {\n",\
"\n",\
"	$hlist add Topology -itemtype text -text \"Topology\"\n",\
"	\n",\
"	set i [expr 0]\n",\
"	while {$i < $nrentities} {\n",\
"	    set entity [lindex $entities [expr $i]]\n",\
"	    incr i 1\n",\
"	    set entityname [lindex $entities [expr $i]]\n",\
"	    $hlist add Topology/$entity -text $entityname -data $entityname\n",\
"	    incr i 1\n",\
"	    $w.mtre close Topology/$entity\n",\
"	}\n",\
"	\n",\
"	$w.mtre autosetmode\n",\
"	$w.mtre open Topology\n",\
"	\n",\
"	set i [expr 0]\n",\
"	while {$i < $nrentities} {\n",\
"	    set entity [lindex $entities [expr $i]]\n",\
"	    $w.mtre close Topology/$entity\n",\
"	    incr i 2\n",\
"	}\n",\
"\n",\
"	$w.mtre autosetmode\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"proc rebuildacisdialog {} {\n",\
"    if {[winfo exists .acis_dlg] == 1} {\n",\
"	[.acis_dlg.mtre subwidget hlist] delete all\n",\
"	acisdialogbuildtree \n",\
"    }\n",\
"}\n",\
"\n",\
"proc checkacisloaded { } {\n",\
"    set isacisgeometryloaded [Ng_ACISCommand isacisgeometryloaded]\n",\
"    if {$isacisgeometryloaded == 0} {\n",\
"	puts \"no IGES/STEP geometry loaded\"\n",\
"	destroy .acis_dlg\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc selectentity { entityname } {\n",\
"    global entities\n",\
"    set nrentities [expr [llength $entities]]\n",\
"    set i [expr 0]\n",\
"    while {$i < $nrentities} {\n",\
"	set entitylength []\n",\
"        \n",\
"	set entity2 [lindex $entities [expr $i]]\n",\
"	incr i 1\n",\
"	set entityname2 [lindex $entities [expr $i]]\n",\
"	incr i 1\n",\
"	set entityname2 [string range $entityname2 0 [expr [string length $entityname]-1]]\n",\
"	\n",\
"	if {$entityname == $entityname2} {\n",\
"	    set hlist [.acis_dlg.mtre subwidget hlist]\n",\
"	    .acis_dlg.mtre open Topology\n",\
"	    set slashpos [string last \"/\" $entity2]\n",\
"	    set entity3 [string range $entity2 0 [expr $slashpos-1]]\n",\
"	    while {$slashpos != -1} {\n",\
"		.acis_dlg.mtre open Topology/$entity3\n",\
"		\n",\
"		set slashpos [string last \"/\" $entity3]\n",\
"		set entity3 [string range $entity3 0 [expr $slashpos-1]]\n",\
"	    }\n",\
"	    $hlist selection clear\n",\
"	    $hlist see Topology/$entity2\n",\
"	    $hlist selection set Topology/$entity2\n",\
"	} \n",\
"    }	    \n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"proc acisdialog { } {\n",\
"    \n",\
"    uplevel 1 {\n",\
"        \n",\
"        global entities\n",\
"        set selectvisual geometry\n",\
"        Ng_SetVisParameters\n",\
"        redraw\n",\
"        \n",\
"        set w .acis_dlg\n",\
"\n",\
"        if {[winfo exists .acis_dlg] == 1} {\n",\
"            wm withdraw $w\n",\
"            wm deiconify $w\n",\
"            focus $w \n",\
"        } {	\n",\
"            toplevel $w\n",\
"            \n",\
"            tixTree $w.mtre -options { separator \"/\" }\n",\
"            pack $w.mtre -fill both -expand yes\n",\
"\n",\
"\n",\
"            acisdialogbuildtree\n",\
"\n",\
"            set hlist [$w.mtre subwidget hlist]\n",\
"\n",\
"\n",\
"            \n",\
"\n",\
"            set solname {\"\"}\n",\
"\n",\
"            puts \"acisdialog2\"\n",\
"            \n",\
"\n",\
"            bind $hlist <ButtonRelease-1> {\n",\
"\n",\
"                set entry [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"                set hentry [string trimleft $entry Topology/]\n",\
"\n",\
"                Ng_ACISCommand selectentity $hentry\n",\
"                redraw\n",\
"            }\n",\
"            \n",\
"\n",\
"            bind $hlist <Double-1> {\n",\
"\n",\
"                puts \"double 1\"\n",\
"\n",\
"                set oldsolname {$solname}\n",\
"                set solname [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"\n",\
"                puts \"solname = $solname\"\n",\
"\n",\
"                if {$solname != \"\" && $oldsolname != $solname } {\n",\
"                    set seppos [string first \"/\" $solname]\n",\
"                    set rootname [string range $solname 0 [expr $seppos-1]]\n",\
"                    \n",\
"                    set entityname [[.acis_dlg.mtre subwidget hlist] info data $solname]\n",\
"                    set spacepos [string first \" \" $entityname]\n",\
"                    set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"                    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"                    set spacepos2 [string first \" \" $helpstring]\n",\
"                    set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"                    if {$rootname == \"Topology\"} {\n",\
"                        Ng_ACISCommand highlightentity $entitytype $entitynumber\n",\
"                        set selectvisual geometry\n",\
"                        redraw\n",\
"                    } {\n",\
"                        set brackpos [string first \" (\" $entityname]\n",\
"                        if {$brackpos != -1} {\n",\
"                            set entityname [string range $entityname 0 $brackpos]\n",\
"                        }\n",\
"\n",\
"                        selectentity $entityname\n",\
"                    }\n",\
"                }\n",\
"            }\n",\
"            \n",\
"            button $w.cl -text \"Close\" -command {\n",\
"                destroy .acis_dlg\n",\
"            }\n",\
"            \n",\
"            puts \"acisdialog3\"\n",\
"\n",\
"            button $w.show -text \"Show\" -command {\n",\
"                set solname [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"                set entityname [[.acis_dlg.mtre subwidget hlist] info data $solname]\n",\
"                set spacepos [string first \" \" $entityname]\n",\
"                set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"                set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"                set spacepos2 [string first \" \" $helpstring]\n",\
"                set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"                Ng_ACISCommand show $entitytype $entitynumber\n",\
"                set selectvisual geometry\n",\
"                                redraw\n",\
"            }\n",\
"            button $w.hide -text \"Hide\" -command {\n",\
"                set solname [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"                set entityname [[.acis_dlg.mtre subwidget hlist] info data $solname]\n",\
"                set spacepos [string first \" \" $entityname]\n",\
"                set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"                set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"                set spacepos2 [string first \" \" $helpstring]\n",\
"                set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"                Ng_ACISCommand hide $entitytype $entitynumber\n",\
"                set selectvisual geometry\n",\
"                                redraw\n",\
"            }\n",\
"\n",\
"            button $w.swaporientation -text \"Swap orientation\" -command {\n",\
"                set solname [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"                set entityname [[.acis_dlg.mtre subwidget hlist] info data $solname]\n",\
"                set spacepos [string first \" \" $entityname]\n",\
"                set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"                set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"                set spacepos2 [string first \" \" $helpstring]\n",\
"                set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"\n",\
"                Ng_ACISCommand swaporientation $entitytype $entitynumber\n",\
"                set selectvisual geometry\n",\
"                                redraw\n",\
"\n",\
"                [.acis_dlg.mtre subwidget hlist] delete all\n",\
"                acisdialogbuildtree	\n",\
"            }\n",\
"\n",\
"            button $w.marksingular -text \"Mark/Unmark as singular\" -command {\n",\
"                set solname [[.acis_dlg.mtre subwidget hlist] info selection]\n",\
"                set entityname [[.acis_dlg.mtre subwidget hlist] info data $solname]\n",\
"                set spacepos [string first \" \" $entityname]\n",\
"                if { $spacepos != 0 } {\n",\
"                    set entitytype [string range $entityname 0 [expr $spacepos-1]]\n",\
"                    set helpstring [string range $entityname [expr $spacepos+1] [expr [string length $entityname]-1]]\n",\
"                    set spacepos2 [string first \" \" $helpstring]\n",\
"                    if { $spacepos2 != 0 } {\n",\
"                        set entitynumber [string range $helpstring 0 [expr $spacepos2-1]]\n",\
"                        \n",\
"                        global ismarkedsingular\n",\
"                        Ng_ACISCommand marksingular $entitytype $entitynumber\n",\
"                        \n",\
"                        set hlist [$w.mtre subwidget hlist]\n",\
"                        \n",\
"                                                                        set style1 [tixDisplayStyle imagetext -foreground black -background white -selectforeground white -selectbackground blue]\n",\
"                        set style2 [tixDisplayStyle imagetext -foreground red -background white -selectforeground red -selectbackground blue]\n",\
"                        \n",\
"                        if { $ismarkedsingular == 0 } {\n",\
"                            $hlist entryconfigure $solname -style $style1\n",\
"                        } {\n",\
"                            $hlist entryconfigure $solname -style $style2\n",\
"                        }\n",\
"\n",\
"                                                                                                                                                                                                                    }\n",\
"                }\n",\
"                \n",\
"                                                                \n",\
"                                            }\n",\
"\n",\
"\n",\
"            checkbutton $w.zoomtohighlightedentity -text \"Zoom to highlighted entity\" \\\n",\
"                -variable acisoptions.zoomtohighlightedentity \\\n",\
"                -command {\n",\
"                    Ng_SetACISVisParameters\n",\
"                    if { ${acisoptions.zoomtohighlightedentity} == 1} {\n",\
"                        set selectvisual geometry\n",\
"                                                Ng_ACISCommand redrawstatus 1\n",\
"                        redraw\n",\
"                    } {\n",\
"                        Ng_ACISCommand redrawstatus 0\n",\
"                    }\n",\
"                }\n",\
"\n",\
"\n",\
"\n",\
"            frame $w.healing -relief groove -borderwidth 3\n",\
"\n",\
"            button $w.healing.checkentities -text \"Analyze geometry\" -command {\n",\
"                set irregent [Ng_ACISCommand findsmallentities]\n",\
"\n",\
"                set w .acis_dlg\n",\
"                set hlist [$w.mtre subwidget hlist]\n",\
"                \n",\
"                $hlist add ProblematicEntities -text \"Problematic Entities\"\n",\
"                $hlist delete offsprings ProblematicEntities\n",\
"\n",\
"                set nritems [expr [llength $irregent]]\n",\
"                set i [expr 0]\n",\
"                while {$i < $nritems} {\n",\
"                    set entity [lindex $irregent [expr $i]]\n",\
"                    incr i 1\n",\
"                    set entityname [lindex $irregent [expr $i]]\n",\
"                    $hlist add ProblematicEntities/$entity -text $entityname -data $entityname\n",\
"                    incr i 1\n",\
"                }\n",\
"                $w.mtre open ProblematicEntities\n",\
"                $w.mtre autosetmode\n",\
"            }\n",\
"\n",\
"            tixControl $w.healing.tolerance -label \"Healing tolerance: \" -integer false \\\n",\
"                -variable acisoptions.tolerance -min 1e-9 -max 1e6 \\\n",\
"                -options {\n",\
"                    entry.width 6\n",\
"                    label.width 25\n",\
"                    label.anchor e\n",\
"                }	\n",\
"\n",\
"            checkbutton $w.healing.fixsmalledges -text \"Fix small edges\" \\\n",\
"                -variable acisoptions.fixsmalledges\n",\
"            \n",\
"            checkbutton $w.healing.fixspotstripfaces -text \"Fix spot/strip faces\" \\\n",\
"                -variable acisoptions.fixspotstripfaces\n",\
"            \n",\
"            checkbutton $w.healing.sewfaces -text \"Sew faces\" \\\n",\
"                -variable acisoptions.sewfaces\n",\
"            \n",\
"            checkbutton $w.healing.makesolids -text \"Make solids\" \\\n",\
"                -variable acisoptions.makesolids\n",\
"            \n",\
"            checkbutton $w.healing.splitpartitions -text \"Split partitions\" \\\n",\
"                -variable acisoptions.splitpartitions\n",\
"            \n",\
"            button $w.healing.heal -text \"Heal geometry\" -command { \n",\
"                .acis_dlg.healing.tolerance invoke\n",\
"                Ng_ACISCommand shapehealing\n",\
"                redraw \n",\
"                [.acis_dlg.mtre subwidget hlist] delete all\n",\
"                acisdialogbuildtree\n",\
"            }\n",\
"\n",\
"            pack $w.healing.checkentities\n",\
"\n",\
"            pack $w.healing.tolerance $w.healing.fixsmalledges \\\n",\
"                $w.healing.fixspotstripfaces $w.healing.sewfaces \\\n",\
"                $w.healing.makesolids $w.healing.splitpartitions -anchor w\n",\
"\n",\
"            pack $w.healing.heal	\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"            pack $w.show $w.hide\n",\
"\n",\
"            pack $w.zoomtohighlightedentity -anchor w\n",\
"                        pack $w.swaporientation\n",\
"            pack $w.marksingular\n",\
"            pack $w.healing -fill x\n",\
"            pack $w.cl\n",\
"            \n",\
"            \n",\
"            wm withdraw $w\n",\
"            wm geom $w +100+100\n",\
"            wm deiconify $w\n",\
"            wm title $w \"IGES/STEP Topology Explorer/Doctor\"\n",\
"            focus .acis_dlg\n",\
"        }\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"if { [catch { NGS_GetData } ] == 0 } { \n",\
"    \n",\
"    set progname \"NGSolve\"\n",\
"    wm title . $progname\n",\
"    \n",\
"    .ngmenu add cascade -label \"Solve\" -menu .ngmenu.solve -underline 1\n",\
"    \n",\
"            \n",\
"    menu .ngmenu.solve\n",\
"    .ngmenu.solve add command -label \"Print Equations\" \\\n",\
"	-command { NGS_PrintRegistered }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    menu .ngmenusolvehelp\n",\
"    .ngmenu.solve add cascade -label \"Help\" -menu .ngmenusolvehelp\n",\
"    \n",\
"    .ngmenusolvehelp add command -label \"Coefficient...\" \\\n",\
"	-command { tk_messageBox -title \"Help\" -message  [ NGS_Help  coefficient ] -type ok }\n",\
"    .ngmenusolvehelp add command -label \"Bilinear-form...\" \\\n",\
"	-command { tk_messageBox -title \"Help\" -message  [ NGS_Help  bilinearform ] -type ok }\n",\
"    .ngmenusolvehelp add command -label \"Linear-form...\" \\\n",\
"	-command { tk_messageBox -title \"Help\" -message  [ NGS_Help  linearform ] -type ok }\n",\
"\n",\
"        .ngmenusolvehelp add cascade -label \"Numprocs...\" -menu .ngmenusolvehelpnp \n",\
"    \n",\
"    .ngmenusolvehelp add command -label \"Latest News...\"  \\\n",\
"	-command { tk_messageBox -title \"Latest News\" -message \\\n",\
"                       { \n",\
"                           06042004 online documentation (JS) \n",\
"                       } -type ok }  ;\n",\
"    \n",\
"\n",\
"\n",\
"    .ngmenu.solve add command -label \"Load PDE...\" -accelerator \"<l><p>\"\\\n",\
"	-command { \n",\
"	    set types { {\"Partial Differential Equation\"   {.pde}	} }\n",\
"	    set file [tk_getOpenFile -filetypes $types]\n",\
"	    if {$file != \"\"} {\n",\
"		AddRecentNGSFile $file;\n",\
"                NGS_LoadPDE  $file;  \n",\
"                set selectvisual mesh;\n",\
"                Ng_SetVisParameters	\n",\
"	    }\n",\
"	}\n",\
"    \n",\
"    .ngmenu.solve add cascade -label \"Recent Files\" -menu .ngmenu.solve.recent \n",\
"    menu .ngmenu.solve.recent\n",\
"\n",\
"    .ngmenu.solve add command -label \"Components...\" \\\n",\
"	-command { componentsdialog }\n",\
"    \n",\
"\n",\
"    .ngmenu.solve add command -label \"Print Report\" \\\n",\
"	-command { NGS_PrintPDE }\n",\
"\n",\
"    .ngmenu.solve add command -label \"Memory Usage\" \\\n",\
"	-command { NGS_PrintMemoryUsage }\n",\
"\n",\
"    .ngmenu.solve add command -label \"Print Timing\" \\\n",\
"	-command { NGS_PrintTiming }\n",\
"\n",\
"\n",\
"\n",\
"	    	    	    	    	    \n",\
"\n",\
"\n",\
"\n",\
"    .ngmenu.solve add command -label \"Solve Recent PDE\" -accelerator \"<s><r>\"\\\n",\
"	-command { \n",\
"	    NGS_LoadPDE  [.ngmenu.solve.recent entrycget 1 -label]\n",\
"	    NGS_SolvePDE\n",\
"	    set selectvisual solution\n",\
"	    Ng_SetVisParameters	\n",\
"\n",\
"\n",\
"	    redraw\n",\
"	}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    button .bubar.pde -text \"Recent\" \\\n",\
"	-command { .ngmenu.solve invoke \"Solve Recent PDE\"; }\n",\
"    pack .bubar.pde -side right\n",\
"\n",\
"    button .bubar.solve -text \"Solve\" \\\n",\
"	-command { .ngmenu.solve invoke \"Solve PDE\"; }\n",\
"    pack .bubar.solve -side right\n",\
"\n",\
"    button .bubar.visualize -text \"Visual\" \\\n",\
"	-command { visual_dialog }\n",\
"    pack .bubar.visualize -side right\n",\
"\n",\
"    .ngmenu.solve add command -label \"Solve PDE\" -accelerator \"<s><p>\"\\\n",\
"	-command {\n",\
"	    \n",\
"	    \n",\
"	    \n",\
"	    NGS_SolvePDE\n",\
"	    set selectvisual solution\n",\
"	    Ng_SetVisParameters	\n",\
"\n",\
"	    Ng_Vis_Set parameters; \n",\
"\n",\
"	    \n",\
"	    redraw\n",\
"	}\n",\
"\n",\
"    .ngmenu.solve add cascade -label \"Solve PDE x\" -menu .ngmenu.solve.solvex\n",\
"    menu .ngmenu.solve.solvex\n",\
"\n",\
"    proc SolveX { num } {\n",\
"	for { set i 1 } { $i <= $num } { incr i } {\n",\
"	    uplevel 1 \"NGS_SolvePDE $i\"\n",\
"	}\n",\
"    }\n",\
"\n",\
"    .ngmenu.solve.solvex  add command -label \"1 Level\" -command { SolveX 1 }\n",\
"    .ngmenu.solve.solvex  add command -label \"2 Level\" -command { SolveX 2 }\n",\
"    .ngmenu.solve.solvex  add command -label \"3 Level\" -command { SolveX 3 }\n",\
"    .ngmenu.solve.solvex  add command -label \"4 Level\" -command { SolveX 4 }\n",\
"    .ngmenu.solve.solvex  add command -label \"5 Level\" -command { SolveX 5 }\n",\
"    .ngmenu.solve.solvex  add command -label \"6 Level\" -command { SolveX 6 }\n",\
"    .ngmenu.solve.solvex  add command -label \"7 Level\" -command { SolveX 7 }\n",\
"    .ngmenu.solve.solvex  add command -label \"8 Level\" -command { SolveX 8 }\n",\
"    .ngmenu.solve.solvex  add command -label \"9 Level\" -command { SolveX 9 }\n",\
"    .ngmenu.solve.solvex  add command -label \"10 Level\" -command { SolveX 10 }\n",\
"    .ngmenu.solve.solvex  add command -label \"11 Level\" -command { SolveX 11 }\n",\
"    .ngmenu.solve.solvex  add command -label \"12 Level\" -command { SolveX 12 }\n",\
"    .ngmenu.solve.solvex  add command -label \"13 Level\" -command { SolveX 13 }\n",\
"    .ngmenu.solve.solvex  add command -label \"14 Level\" -command { SolveX 14 }\n",\
"    .ngmenu.solve.solvex  add command -label \"15 Level\" -command { SolveX 15 }\n",\
"    .ngmenu.solve.solvex  add command -label \"16 Level\" -command { SolveX 16 }\n",\
"    .ngmenu.solve.solvex  add command -label \"17 Level\" -command { SolveX 17 }\n",\
"    .ngmenu.solve.solvex  add command -label \"18 Level\" -command { SolveX 18 }\n",\
"    .ngmenu.solve.solvex  add command -label \"19 Level\" -command { SolveX 19 }\n",\
"    .ngmenu.solve.solvex  add command -label \"20 Level\" -command { SolveX 20 }\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"                        \n",\
"\n",\
"    .ngmenu.solve add command -label \"Visualization...\" \\\n",\
"	-command { \n",\
"	    visual_dialog;\n",\
"	}\n",\
"    \n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    \n",\
"    .ngmenu.solve add command -label \"Save Solution...\" \\\n",\
"	-command { \n",\
"	    set types { {\"Solution File\"  {.sol} } }\n",\
"	    set file [tk_getSaveFile -filetypes $types -defaultextension \".sol\"  ]\n",\
"	    if {$file != \"\"} {\n",\
"		NGS_SaveSolution $file \n",\
"	    }\n",\
"	}\n",\
"    \n",\
"    .ngmenu.solve add command -label \"Load Solution...\" \\\n",\
"	-command { \n",\
"	    set types { {\"Solution File\"  {.sol} } }\n",\
"	    set file [tk_getOpenFile -filetypes $types -defaultextension \".sol\"  ]\n",\
"	    if {$file != \"\"} {\n",\
"		NGS_LoadSolution $file \n",\
"		set selectvisual solution\n",\
"		Ng_SetVisParameters\n",\
"		redraw\n",\
"	    }\n",\
"	}\n",\
"\n",\
"    \n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
".ngmenu.help delete \"About...\"\n",\
".ngmenu.help add command -label \"About...\" \\\n",\
"    -command {\n",\
"        tk_messageBox -message \\\n",\
"            \"This is NETGEN/NGSolve \\n mainly written by \\n Joachim Schöberl \\n\\\n",\
"                 at RWTH Aachen University, Germany \\n\\\n",\
"                 and Johannes Kepler University, Linz, Austria \\n\\\n",\
"                 supported by the Austrian Science Foundation FWF \\n\\\n",\
"                 thanks to \\n\\\n",\
"                 F. Bachinger, A. Becirovic, H. Egger, R. Gaisbauer, J. Gerstmayr, U. Langer, A. Sinwel, M. Wabro, S. Zaglmayr\"\n",\
"	}\n",\
"    \n",\
"\n",\
"\n",\
"    proc AddRecentNGSFile { filename } {\n",\
"	global progname\n",\
"	catch { [.ngmenu.solve.recent delete $filename] }\n",\
"	.ngmenu.solve.recent insert 0 command -label $filename \\\n",\
"	    -command \"AddRecentNGSFile {$filename}; \n",\
"		NGS_LoadPDE  {$filename};\n",\
"                set selectvisual mesh;\n",\
"                Ng_SetVisParameters	\n",\
"                wm title . [concat \\\" $progname - $filename \\\"];\"\n",\
"	\n",\
"	if { [.ngmenu.solve.recent index last] >= 6 } {\n",\
"	    .ngmenu.solve.recent delete last }\n",\
"	\n",\
"	savengsinifile;\n",\
"    }\n",\
"    \n",\
"\n",\
"        proc savengsinifile { } {\n",\
"	uplevel 1  {\n",\
"	    set datei [open ngs.ini w]\n",\
"	    for { set i [.ngmenu.solve.recent index last] } { $i >= 1 } { incr i -1 } {\n",\
"		puts $datei \"recentfile \\\"[.ngmenu.solve.recent entrycget $i -label]\\\"\"\n",\
"	    }\n",\
"	\n",\
"	    close $datei\n",\
"	}\n",\
"    }\n",\
"    \n",\
"    proc loadngsinifile { } {\n",\
"	if { [file exists ngs.ini] == 1 } {\n",\
"	    set datei [open ngs.ini r]\n",\
"	    while { [gets $datei line] >= 0 } {\n",\
"		if {[lindex $line 0] == \"recentfile\"} {\n",\
"		    AddRecentNGSFile [lindex $line 1]\n",\
"		}\n",\
"	    }\n",\
"	    close $datei\n",\
"	}\n",\
"    }\n",\
"\n",\
"\n",\
"loadngsinifile;\n",\
"    \n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"    proc componentsdialog { } {\n",\
"	\n",\
"	set w .components_dlg\n",\
"	\n",\
"	if {[winfo exists .components_dlg] == 1} {\n",\
"	    wm withdraw $w\n",\
"	    wm deiconify $w\n",\
"	    focus $w \n",\
"	} {\n",\
"\n",\
"	    toplevel $w\n",\
"\n",\
"\n",\
"	    tixTree $w.mtre -options { separator \"\\\\\" }\n",\
"	    pack $w.mtre -fill both -expand y\n",\
"	    set hlist [$w.mtre subwidget hlist]\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	    $hlist add coeffs -itemtype text -text \"Coefficients\"\n",\
"	    set coefs [NGS_GetData coefficients]\n",\
"	    foreach coef $coefs {\n",\
"		$hlist add coeffs\\\\$coef -itemtype text -text $coef\n",\
"	    }\n",\
"\n",\
"\n",\
" 	    $hlist add spaces -itemtype text -text \"Spaces\"\n",\
"	    set spaces [NGS_GetData spaces]\n",\
"	    foreach space $spaces {\n",\
"		$hlist add spaces\\\\$space -itemtype text -text $space\n",\
"	    }\n",\
"\n",\
" 	    $hlist add biforms -itemtype text -text \"Bilinear-forms\"\n",\
"	    set biforms [NGS_GetData bilinearforms]\n",\
"	    foreach biform $biforms {\n",\
"		$hlist add biforms\\\\$biform -itemtype text -text $biform\n",\
"	    }\n",\
"\n",\
" 	    $hlist add liforms -itemtype text -text \"Linear-forms\"\n",\
"	    set liforms [NGS_GetData linearforms]\n",\
"	    foreach liform $liforms {\n",\
"		$hlist add liforms\\\\$liform -itemtype text -text $liform\n",\
"	    }\n",\
"\n",\
" 	    $hlist add gridfuns -itemtype text -text \"Grid-functions\"\n",\
"	    set gridfuns [NGS_GetData gridfunctions]\n",\
"	    foreach gridfun $gridfuns {\n",\
"		$hlist add gridfuns\\\\$gridfun -itemtype text -text $gridfun\n",\
"	    }\n",\
"\n",\
" 	    $hlist add preconds -itemtype text -text \"Preconditioners\"\n",\
"	    set preconds [NGS_GetData preconditioners]\n",\
"	    foreach precond $preconds {\n",\
"		$hlist add preconds\\\\$precond -itemtype text -text $precond\n",\
"	    }\n",\
"\n",\
" 	    $hlist add numprocs -itemtype text -text \"NumProcs\"\n",\
"	    set numprocs [NGS_GetData numprocs]\n",\
"	    foreach numproc $numprocs {\n",\
"		$hlist add numprocs\\\\$numproc -itemtype text -text $numproc\n",\
"	    }\n",\
"\n",\
"\n",\
"\n",\
"	    \n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"	    $w.mtre autosetmode\n",\
"	    \n",\
"\n",\
"	    bind $hlist <Double-1> {\n",\
"		set solname [[.components_dlg.mtre subwidget hlist] info selection]\n",\
"		puts $solname\n",\
"		set seppos [string first \\\\ $solname]\n",\
"		if { $seppos != -1 } {\n",\
"		    set field [string range $solname 1 [expr $seppos-1]]\n",\
"		    set name [string range $solname [expr $seppos+1] [expr [string length $solname]-2]]\n",\
"		    puts \"field = $field, name = $name\"\n",\
"		    NGS_PrintPDE $field $name\n",\
"		}\n",\
"	    }\n",\
"\n",\
"	    button $w.cl -text \"Close\" -command {\n",\
"		destroy .components_dlg\n",\
"	    }\n",\
"\n",\
"	    pack  $w.cl\n",\
"	    \n",\
"	    \n",\
"	    wm withdraw $w\n",\
"	    wm geom $w +100+100\n",\
"	    wm deiconify $w\n",\
"	    wm title $w \"Components\"\n",\
"	    focus .components_dlg\n",\
"	}\n",\
"    }\n",\
"\n",\
"bind . <l><p> { .ngmenu.solve invoke \"Load PDE...\" }  ; \n",\
"bind . <s><r> { .ngmenu.solve invoke \"Solve Recent PDE\" }  ; \n",\
"bind . <s><p> { .ngmenu.solve invoke \"Solve PDE\" }  ; \n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"set zugstange 0\n",\
"\n",\
"\n",\
"catch { source ${ngdir}/trafo/menu.tcl }\n",\
"\n",\
"\n",\
"\n",\
"setgranularity ${meshoptions.fineness}\n",\
"\n",\
"Ng_SetMeshingParameters\n",\
"Ng_SetVisParameters\n",\
"Ng_SetDebugParameters\n",\
"Ng_STLDoctor\n",\
"Ng_GeometryOptions set\n",\
"Ng_SetOCCVisParameters\n",\
"\n",\
"if { $batchmode != \"defined\" } {\n",\
"    catch { \n",\
"	wm protocol . WM_DELETE_WINDOW { .ngmenu.file invoke \"Quit\" }\n",\
"	wm deiconify .\n",\
"    }\n",\
"}\n",\
"\n",\
"set trafoapp 0\n",\
"catch { source ${ngdir}/trafoapp/trafoapp.tcl }\n",\
"\n",\
"set geofilename [Ng_GetCommandLineParameter geofile]\n",\
"\n",\
"if { $geofilename != \"undefined\" && \n",\
"     [info exists trafo] == 0 && $zugstange == 0} {\n",\
"\n",\
"    if { [ catch { Ng_LoadGeometry $geofilename } errstring] == 0 } {\n",\
"	if { $batchmode != \"defined\" } {\n",\
"	    AddRecentFile $geofilename\n",\
"	}\n",\
"	Ng_ParseGeometry\n",\
"	if { $batchmode != \"defined\" } {\n",\
"	    set selectvisual geometry\n",\
"	    Ng_SetVisParameters\n",\
"	    redraw\n",\
"	    wm title . [concat \"$progname - \" $geofilename]\n",\
"	}\n",\
"	set dirname [file dirname $geofilename]\n",\
"	set basefilename [file tail [file rootname $geofilename]]\n",\
"    } {\n",\
"	puts \"Problem with input file:\"\n",\
"	puts \"$errstring\"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"set cnt 0\n",\
"foreach { gran } { verycoarse coarse moderate fine veryfine } {\n",\
"    set cnt [expr $cnt + 1]\n",\
"    if { [Ng_GetCommandLineParameter $gran] == \"defined\" } {\n",\
"	set meshoptions.fineness $cnt\n",\
"	setgranularity ${meshoptions.fineness}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"set meshfilename [Ng_GetCommandLineParameter meshfile]\n",\
"if { $meshfilename == \"undefined\" } {\n",\
"    set meshfilename out.mesh\n",\
"}\n",\
"\n",\
"set meshfiletype [Ng_GetCommandLineParameter meshfiletype]\n",\
"if { $meshfiletype == \"undefined\" } {\n",\
"    set meshfiletype netgen\n",\
"}\n",\
"\n",\
"set inputmeshfilename [Ng_GetCommandLineParameter inputmeshfile]\n",\
"\n",\
"set mergemeshfilename [Ng_GetCommandLineParameter mergefile]\n",\
"\n",\
"set meshsizefilename [Ng_GetCommandLineParameter meshsizefile]\n",\
"\n",\
"if { $meshsizefilename != \"undefined\" } {\n",\
"    set options.meshsizefilename $meshsizefilename\n",\
"}\n",\
"\n",\
"set refinementfilename [Ng_GetCommandLineParameter refinementfile]\n",\
"\n",\
"\n",\
"if { $batchmode == \"defined\" && $solvemode != \"defined\"} {\n",\
"    set options.parthread 0\n",\
"    if { $shellmode == \"undefined\" } {\n",\
"      set selectvisual mesh\n",\
"      Ng_SetVisParameters\n",\
"\n",\
"      set meshsize [Ng_GetCommandLineParameter meshsize]\n",\
"      if {$meshsize != \"undefined\"} { set options.meshsize $meshsize }\n",\
"        \n",\
"      if { $inputmeshfilename == \"undefined\" } {\n",\
"	Ng_GenerateMesh ${meshoptions.firststep} ${meshoptions.laststep}\n",\
"      } else {\n",\
"	Ng_LoadMesh $inputmeshfilename\n",\
"	if { $mergemeshfilename != \"undefined\" } {\n",\
"	    Ng_MergeMesh $mergemeshfilename\n",\
"        }\n",\
"      }\n",\
"	\n",\
"      if { $refinementfilename != \"undefined\" } {\n",\
"	  Ng_Bisect $refinementfilename\n",\
"      }\n",\
"\n",\
"      if { $meshfiletype == \"netgen\" } {\n",\
"	Ng_SaveMesh $meshfilename\n",\
"      } else {\n",\
"	if { [catch { Ng_ExportMesh $meshfilename $meshfiletype } ] == 1 } {\n",\
"	    puts \"Unknown file format $meshfiletype\"\n",\
"        }\n",\
"      }\n",\
"      Ng_Exit;\n",\
"\n",\
"      exit\n",\
"  } else {\n",\
"      set code [catch { source ${ngdir}/ngtcltk/ngshell.tcl } errcode]\n",\
"      if {$code} {\n",\
"	  puts \"error: $errcode\"\n",\
"      }  \n",\
"      set code [ catch {Ng_RunShell} errcode]\n",\
"      if {$code} {\n",\
"	  puts \"error: $errcode\"\n",\
"      }  \n",\
"      \n",\
"      Ng_Exit;\n",\
"      exit\n",\
"  }\n",\
"    \n",\
"}\n",\
"\n",\
"set stereo [Ng_GetCommandLineParameter stereo]\n",\
"if { $stereo == \"defined\" } {\n",\
"    set viewoptions.stereo 1 \n",\
"    puts \"use stereo mode\" \n",\
"    Ng_SetVisParameters; \n",\
"    redraw \n",\
"}\n",\
"\n",\
"\n",\
"set scriptfilename [Ng_GetCommandLineParameter script]\n",\
"if { $scriptfilename != \"undefined\" } {\n",\
"    if { [catch { source $scriptfilename } errstring] == 1 } {\n",\
"	puts \"Error in input: $errstring\"\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"if { [Ng_GetCommandLineParameter help]==\"defined\" } {\n",\
"    if { $zugstange == 1 } {\n",\
"	print_zug_commandline_help\n",\
"	exit;\n",\
"    } {\n",\
"	if { $trafoapp == 1 } {\n",\
"	    print_trafo_commandline_help;\n",\
"	} {\n",\
"	    print_commandline_help; \n",\
"	    Ng_Exit;\n",\
"	    exit\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"if { [file exists startup.tcl] } {\n",\
"    source startup.tcl }\n",\
"\n",\
"if { [Ng_GetCommandLineParameter recent]==\"defined\" } {\n",\
"    if { [catch { .ngmenu.solve invoke \"Solve Recent PDE\";  } errstring] == 1 } {\n",\
"	puts \"TCL-ERROR handler:\\n $errstring\";\n",\
"	exit;\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"set pdefilename [Ng_GetCommandLineParameter pdefile]\n",\
"if { $pdefilename != \"undefined\" } {\n",\
"    NGS_LoadPDE  $pdefilename;  \n",\
"\n",\
"    set solve [Ng_GetCommandLineParameter solve]\n",\
"    if { $zugstange == 1 } {\n",\
"	set options.parthread 0\n",\
"	NGS_SolvePDE;\n",\
"    } {\n",\
"	if { $solve == \"defined\" } {\n",\
"	    set options.parthread 0\n",\
"	    NGS_SolvePDE\n",\
"	    exit;\n",\
"	} {\n",\
"	    if { $solve != \"undefined\" } {\n",\
"		set options.parthread 0\n",\
"		for { set l 1 } { $l <= $solve } { incr l } { NGS_SolvePDE $l }\n",\
"		exit;\n",\
"	    }\n",\
"	}\n",\
"    }\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"catch { source ${ngdir}/trafo/trafo.tcl }\n",\
"\n",\
"catch { source ${ngdir}/trafoapp/smallmodels.tcl }\n",\
"\n",\
"catch { \n",\
"  source ${ngdir}/ngtcltk/ngshell.tcl\n",\
"  source ${ngdir}/ngtcltk/ngtesting.tcl\n",\
"}\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
0};
