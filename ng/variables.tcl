# netgen global tcl-variables

set drawmode rotate
set selectvisual geometry

set dirname .
set loadgeomtypevar "All Geometry types"

set basefilename filename

set meshoptions.fineness 3
set meshoptions.firststep ag
set meshoptions.laststep ov
set options.memory 0

set options.localh 1
set options.delaunay 1
set options.checkoverlap 1
set options.checkoverlappingboundary 0
set options.checkchartboundary 1
set options.startinsurface 0
set options.blockfill 1
set options.debugmode 0
set options.dooptimize 1
set options.parthread 1
set options.elsizeweight 0.2
set options.secondorder 0
set options.elementorder 1
set options.quad 0
set options.inverttets 0
set options.inverttrigs 0
set options.autozrefine 0


set options.meshsize 1000
set options.minmeshsize 0

set options.curvaturesafety 2
set options.segmentsperedge 2
set options.meshsizefilename ""
set options.badellimit 175
set options.optsteps2d 3
set options.optsteps3d 5
set options.opterrpow 2

set options.grading 0.5
set options.printmsg 2

set debug.slowchecks 0
set debug.debugoutput 0
set debug.haltexistingline 0
set debug.haltoverlap 0
set debug.haltsuccess 0
set debug.haltnosuccess 0
set debug.haltlargequalclass 0
set debug.haltsegment 0
set debug.haltnode 0
set debug.haltface 0
set debug.haltfacenr 0
set debug.haltsegmentp1 0
set debug.haltsegmentp2 0

set geooptions.drawcsg 1
set geooptions.detail 0.001
set geooptions.accuracy 1e-6
set geooptions.facets 20
set geooptions.minx -1000
set geooptions.miny -1000
set geooptions.minz -1000
set geooptions.maxx 1000
set geooptions.maxy 1000
set geooptions.maxz 1000

set viewqualityplot 0
set memuseplot 0
set viewrotatebutton 0
set showsensitivehelp 0
set showhelpline 0

set viewoptions.specpointvlen 0.3
set viewoptions.light.amb 0.3
set viewoptions.light.diff 0.7
set viewoptions.light.spec 1
set viewoptions.light.locviewer 0
set viewoptions.mat.shininess 50
set viewoptions.mat.transp 0.3
set viewoptions.colormeshsize 0
set viewoptions.whitebackground 1
set viewoptions.drawcoordinatecross 1
set viewoptions.drawcolorbar 1
set viewoptions.drawnetgenlogo 1
set viewoptions.stereo 0
set viewoptions.shrink 1

set viewoptions.drawfilledtrigs 1
set viewoptions.drawedges 0
set viewoptions.drawbadels 0
set viewoptions.centerpoint 0
set viewoptions.drawelement 0
set viewoptions.drawoutline 1
set viewoptions.drawtets 0
set viewoptions.drawtetsdomain 0
set viewoptions.drawprisms 0
set viewoptions.drawpyramids 0
set viewoptions.drawhexes 0
set viewoptions.drawidentified 0
set viewoptions.drawpointnumbers 0
set viewoptions.drawedgenumbers 0
set viewoptions.drawfacenumbers 0
set viewoptions.drawelementnumbers 0
set viewoptions.drawdomainsurf 0

set viewoptions.drawededges 1
set viewoptions.drawedpoints 1
set viewoptions.drawedpointnrs 0
set viewoptions.drawedtangents 0
set viewoptions.drawededgenrs 0
set viewoptions.drawmetispartition 0

set viewoptions.drawcurveproj 0
set viewoptions.drawcurveprojedge 1

set viewoptions.clipping.nx 0
set viewoptions.clipping.ny 1
set viewoptions.clipping.nz 0
set viewoptions.clipping.dist 0
set viewoptions.clipping.dist2 0
set viewoptions.clipping.enable 0
set viewoptions.clipping.onlydomain 0
set viewoptions.clipping.notdomain 0

set viewoptions.usecentercoords 0
set viewoptions.centerx 0
set viewoptions.centery 0
set viewoptions.centerz 0

set viewoptions.drawspecpoint 0
set viewoptions.specpointx 0
set viewoptions.specpointy 0
set viewoptions.specpointz 0


set stloptions.showtrias 0
set stloptions.showfilledtrias 1
set stloptions.showedges 1
set stloptions.showmarktrias 0
set stloptions.showactivechart 0
set stloptions.yangle 30
set stloptions.contyangle 20
set stloptions.edgecornerangle 60
set stloptions.chartangle 15
set stloptions.outerchartangle 70
set stloptions.usesearchtree 0
set stloptions.chartnumber 1
set stloptions.charttrignumber 1
set stloptions.chartnumberoffset 0

set stloptions.atlasminh 0.1
set stloptions.resthsurfcurvfac 2
set stloptions.resthsurfcurvenable 0
set stloptions.resthatlasfac 2
set stloptions.resthatlasenable 1
set stloptions.resthchartdistfac 1.2
set stloptions.resthchartdistenable 1
set stloptions.resthlinelengthfac 0.5
set stloptions.resthlinelengthenable 1
set stloptions.resthcloseedgefac 1
set stloptions.resthcloseedgeenable 1
set stloptions.resthminedgelen 0.01
set stloptions.resthminedgelenenable 1
set stloptions.resthedgeanglefac 1
set stloptions.resthedgeangleenable 0
set stloptions.resthsurfmeshcurvfac 1
set stloptions.resthsurfmeshcurvenable 0
set stloptions.recalchopt 1

set stldoctor.drawmeshededges 1
set stldoctor.geom_tol_fact 0.000001
set stldoctor.useexternaledges 0
set stldoctor.showfaces 0
set stldoctor.conecheck 1
set stldoctor.spiralcheck 1
set stldoctor.selecttrig 0
set stldoctor.selectmode 1
set stldoctor.longlinefact 0
set stldoctor.showexcluded 1
set stldoctor.edgeselectmode 0
set stldoctor.nodeofseltrig 1
set stldoctor.showtouchedtrigchart 0
set stldoctor.showedgecornerpoints 0
set stldoctor.showmarkedtrigs 1
set stldoctor.dirtytrigfact 0.01
set stldoctor.smoothangle 90
set stldoctor.selectwithmouse 1
set stldoctor.showvicinity 0
set stldoctor.vicinity 50
set stldoctor.smoothnormalsweight 0.2

set occoptions.showvolumenr 0
set occoptions.showsurfaces 1
set occoptions.showedges 1
set occoptions.showsolidnr 0
set occoptions.showsolidnr2 0
set occoptions.visproblemfaces 0
set occoptions.zoomtohighlightedentity 0
set occoptions.deflection 1
set occoptions.tolerance 1e-3
set occoptions.fixsmalledges 1
set occoptions.fixspotstripfaces 1
set occoptions.sewfaces 1
set occoptions.makesolids 1
set occoptions.splitpartitions 0

set meshdoctor.active 0
set meshdoctor.markedgedist 1


# variablenname mit punkt problematisch!
set status_np 0
set status_ne 0
set status_nse 0
set status_working " "
set status_task " "
set status_percent 0
set status_filename 0
set status_tetqualclasses "10 20 30 40 10 20 30 40 10 20 30 40 10 20 30 40 10 20 30 40"

set exportfiletype "Neutral Format"

set preproc.facenr 0
set preproc.selectmode query
set preproc.numtrig 0

set mem_moveable 0


set multithread_pause 0
set multithread_testmode 0
set multithread_redraw 0
set multithread_drawing 0
set multithread_terminate 0
set multithread_running 0

set level 0


set tablesforoutput {}



set optlist {
    options.localh 
    options.delaunay 
    options.checkoverlap 
    options.startinsurface 
    options.blockfill 
    options.dooptimize 
    options.elsizeweight 
    options.meshsize 
    options.minmeshsize 
    options.curvaturesafety 
    options.optsteps2d 
    options.optsteps3d 
    options.secondorder
}


set visoptions.usetexture 1
set visoptions.invcolor 0
set visoptions.imaginary 0
set visoptions.lineartexture 0
set visoptions.numtexturecols 16
set visoptions.showclipsolution 1
set visoptions.showsurfacesolution 0
set visoptions.drawfieldlines 0
set visoptions.drawpointcurves 1
set visoptions.numfieldlines 100
set visoptions.fieldlinesrandomstart 0
set visoptions.fieldlinesstartarea box
set visoptions.fieldlinesstartareap1x 1
set visoptions.fieldlinesstartareap1y 1
set visoptions.fieldlinesstartareap1z 1
set visoptions.fieldlinesstartareap2x 0
set visoptions.fieldlinesstartareap2y 0
set visoptions.fieldlinesstartareap2z 0
set visoptions.fieldlinesstartface -1
set visoptions.fieldlinesfilename none
set visoptions.fieldlinestolerance 0.0005
set visoptions.fieldlinesrktype crungekutta
set visoptions.fieldlineslength 0.5
set visoptions.fieldlinesmaxpoints 500
set visoptions.fieldlinesthickness 0.0015
set visoptions.fieldlinesvecfunction none
set visoptions.fieldlinesphase 0
set visoptions.fieldlinesonlyonephase 1


set visoptions.lineplotfile empty
set visoptions.lineplotsource file
set visoptions.lineplotusingx 0
set visoptions.lineplotusingy 1
set visoptions.lineplotautoscale 1
set visoptions.lineplotxmin 0
set visoptions.lineplotxmax 1
set visoptions.lineplotymin 0
set visoptions.lineplotymax 1
set visoptions.lineplotcurrentnum -1
set visoptions.lineplotinfos ""
set visoptions.lineplotselected none
set visoptions.lineplotselector ""
set visoptions.lineplotcolor red
set visoptions.lineplotsizex 500
set visoptions.lineplotsizey 400
set visoptions.lineplotselectedeval 0
set visoptions.lineplotdatadescr "column1 column2 column3"
set visoptions.lineplotxcoordselector ""
set visoptions.lineplotycoordselector ""
set visoptions.evaluatefilenames none
set visoptions.evaluatefiledescriptions none


set visoptions.clipsolution none
set visoptions.scalfunction none
set visoptions.vecfunction none
set visoptions.evaluate abs
set visoptions.gridsize 20
set visoptions.xoffset 0
set visoptions.yoffset 0
set visoptions.autoscale 1
set visoptions.redrawperiodic 0
set visoptions.logscale 0
set visoptions.mminval 0
set visoptions.mmaxval 1
set visoptions.isolines 0
set visoptions.isosurf 0
set visoptions.subdivisions 1
set visoptions.numiso 10
set visoptions.autoredraw 0
set visoptions.autoredrawtime 2
set visoptions.simulationtime 0
set visoptions.multidimcomponent 0

# deform by vector function
set visoptions.deformation 0
set visoptions.scaledeform1 1
set visoptions.scaledeform2 1

set parallel_netgen 0















set optfilename [file join $nguserdir ng.opt]
set inifilename [file join $nguserdir ng.ini]
set meshinifilename [file join $nguserdir ngmesh.ini]

global env
if { [llength [array names env NG_OPT]] == 1 } {
    if { [string length $env(NG_OPT)] > 0 } {
	set optfilename $env(NG_OPT) 
    }
}

if { [file exists $optfilename] == 1 } {
    set datei [open $optfilename r]
    while { [gets $datei line] >= 0 } {
	set [lindex $line 0] [lindex $line 1]
    }
    close $datei
} {
    puts "optfile $optfilename does not exist - using default values"
}




proc saveoptions { } {
    uplevel 1  {
	set file $optfilename
	
	if {$file != ""} {
	    set datei [open $file w]
	    puts $datei "dirname  ${dirname}"
	    puts $datei "loadgeomtypevar  \"${loadgeomtypevar}\""
	    puts $datei "exportfiletype  \"${exportfiletype}\""
	    puts $datei "meshoptions.fineness  ${meshoptions.fineness}"
	    puts $datei "meshoptions.firststep ${meshoptions.firststep}"
	    puts $datei "meshoptions.laststep  ${meshoptions.laststep}" 
	    puts $datei "options.localh  ${options.localh}"
	    puts $datei "options.delaunay  ${options.delaunay}"
	    puts $datei "options.checkoverlap  ${options.checkoverlap}"
	    puts $datei "options.checkchartboundary  ${options.checkchartboundary}"
	    puts $datei "options.startinsurface  ${options.startinsurface}" 
	    puts $datei "options.blockfill  ${options.blockfill}" 
	    puts $datei "options.debugmode  ${options.debugmode}" 
	    puts $datei "options.dooptimize ${options.dooptimize}" 
	    puts $datei "options.parthread  ${options.parthread}"  
	    puts $datei "options.elsizeweight  ${options.elsizeweight}" 
	    puts $datei "options.secondorder  ${options.secondorder}" 
	    puts $datei "options.elementorder  ${options.elementorder}" 
#	    puts $datei "options.memory  ${options.memory}" 
	    puts $datei "options.quad  ${options.quad}" 
	    puts $datei "options.inverttets  ${options.inverttets}" 
	    puts $datei "options.inverttrigs  ${options.inverttrigs}" 
	    puts $datei "options.autozrefine ${options.autozrefine}" 
	    puts $datei "options.meshsize  ${options.meshsize}" 
	    puts $datei "options.minmeshsize  ${options.minmeshsize}" 
	    puts $datei "options.curvaturesafety  ${options.curvaturesafety}" 
	    puts $datei "options.segmentsperedge  ${options.segmentsperedge}" 
	    puts $datei "options.meshsizefilename  ${options.meshsizefilename}" 
	    puts $datei "options.badellimit  ${options.badellimit}" 
	    puts $datei "options.optsteps2d  ${options.optsteps2d}" 
	    puts $datei "options.optsteps3d  ${options.optsteps3d}" 
	    puts $datei "options.opterrpow  ${options.opterrpow}" 
	    puts $datei "options.grading  ${options.grading}" 
	    puts $datei "options.printmsg  ${options.printmsg}" 
	    puts $datei "geooptions.drawcsg  ${geooptions.drawcsg}" 
	    puts $datei "geooptions.detail  ${geooptions.detail}" 
	    puts $datei "geooptions.accuracy  ${geooptions.accuracy}" 
	    puts $datei "geooptions.facets  ${geooptions.facets}" 
	    puts $datei "geooptions.minx  ${geooptions.minx}" 
	    puts $datei "geooptions.miny  ${geooptions.miny}" 
	    puts $datei "geooptions.minz  ${geooptions.minz}" 
	    puts $datei "geooptions.maxx  ${geooptions.maxx}" 
	    puts $datei "geooptions.maxy  ${geooptions.maxy}" 
	    puts $datei "geooptions.maxz  ${geooptions.maxz}" 
	    puts $datei "viewoptions.specpointvlen  ${viewoptions.specpointvlen}" 
	    puts $datei "viewoptions.light.amb  ${viewoptions.light.amb}" 
	    puts $datei "viewoptions.light.diff ${viewoptions.light.diff}"
	    puts $datei "viewoptions.light.spec ${viewoptions.light.spec}"
	    puts $datei "viewoptions.light.locviewer ${viewoptions.light.locviewer}"
	    puts $datei "viewoptions.mat.shininess  ${viewoptions.mat.shininess}" 
	    puts $datei "viewoptions.mat.transp  ${viewoptions.mat.transp}" 
	    puts $datei "viewoptions.colormeshsize ${viewoptions.colormeshsize}"
	    puts $datei "viewoptions.whitebackground  ${viewoptions.whitebackground}" 
	    puts $datei "viewoptions.drawcolorbar  ${viewoptions.drawcolorbar}" 
	    puts $datei "viewoptions.drawcoordinatecross  ${viewoptions.drawcoordinatecross}" 
	    puts $datei "viewoptions.drawnetgenlogo  ${viewoptions.drawnetgenlogo}" 
	    puts $datei "viewoptions.stereo  ${viewoptions.stereo}" 
	    puts $datei "viewoptions.drawfilledtrigs  ${viewoptions.drawfilledtrigs}" 
	    puts $datei "viewoptions.drawedges  ${viewoptions.drawedges}" 
	    puts $datei "viewoptions.drawbadels  ${viewoptions.drawbadels}" 
	    puts $datei "viewoptions.centerpoint  ${viewoptions.centerpoint}" 
	    puts $datei "viewoptions.drawelement  ${viewoptions.drawelement}" 
	    puts $datei "viewoptions.drawoutline  ${viewoptions.drawoutline}" 
	    puts $datei "viewoptions.drawtets  ${viewoptions.drawtets}"
	    puts $datei "viewoptions.drawprisms  ${viewoptions.drawprisms}"
	    puts $datei "viewoptions.drawpyramids  ${viewoptions.drawpyramids}" 
	    puts $datei "viewoptions.drawhexes  ${viewoptions.drawhexes}" 
	    puts $datei "viewoptions.drawidentified  ${viewoptions.drawidentified}" 
	    puts $datei "viewoptions.drawpointnumbers  ${viewoptions.drawpointnumbers}" 
	    
	    puts $datei "viewoptions.drawededges  ${viewoptions.drawededges}" 
	    puts $datei "viewoptions.drawedpoints  ${viewoptions.drawedpoints}" 
	    puts $datei "viewoptions.drawedpointnrs  ${viewoptions.drawedpointnrs}" 
	    puts $datei "viewoptions.drawedtangents  ${viewoptions.drawedtangents}" 
	    puts $datei "viewoptions.shrink  ${viewoptions.shrink}" 
	    
	    puts $datei "stloptions.showtrias  ${stloptions.showtrias}" 
	    puts $datei "stloptions.showfilledtrias  ${stloptions.showfilledtrias}" 
	    puts $datei "stloptions.showedges  ${stloptions.showedges}" 
	    puts $datei "stloptions.showmarktrias  ${stloptions.showmarktrias}" 
	    puts $datei "stloptions.showactivechart  ${stloptions.showactivechart}" 
	    puts $datei "stloptions.yangle  ${stloptions.yangle}" 
	    puts $datei "stloptions.contyangle  ${stloptions.contyangle}" 
	    puts $datei "stloptions.edgecornerangle  ${stloptions.edgecornerangle}" 
	    puts $datei "stloptions.chartangle  ${stloptions.chartangle}" 
	    puts $datei "stloptions.outerchartangle  ${stloptions.outerchartangle}" 
	    puts $datei "stloptions.usesearchtree  ${stloptions.usesearchtree}" 
	    puts $datei "stloptions.chartnumber  ${stloptions.chartnumber}" 
	    puts $datei "stloptions.charttrignumber  ${stloptions.charttrignumber}" 
	    puts $datei "stloptions.chartnumberoffset  ${stloptions.chartnumberoffset}" 
	    puts $datei "stloptions.atlasminh  ${stloptions.atlasminh}" 
	    puts $datei "stloptions.resthsurfcurvfac  ${stloptions.resthsurfcurvfac}" 
	    puts $datei "stloptions.resthsurfcurvenable  ${stloptions.resthsurfcurvenable}" 
	    puts $datei "stloptions.resthatlasfac  ${stloptions.resthatlasfac}" 
	    puts $datei "stloptions.resthatlasenable  ${stloptions.resthatlasenable}" 
	    puts $datei "stloptions.resthchartdistfac  ${stloptions.resthchartdistfac}" 
	    puts $datei "stloptions.resthchartdistenable  ${stloptions.resthchartdistenable}" 
	    puts $datei "stloptions.resthlinelengthfac  ${stloptions.resthlinelengthfac}" 
	    puts $datei "stloptions.resthlinelengthenable  ${stloptions.resthlinelengthenable}" 
		puts $datei "stloptions.resthminedgelen ${stloptions.resthminedgelen}"
		puts $datei "stloptions.resthminedgelenenable ${stloptions.resthminedgelenenable}"
	    puts $datei "stloptions.resthcloseedgefac  ${stloptions.resthcloseedgefac}" 
	    puts $datei "stloptions.resthcloseedgeenable  ${stloptions.resthcloseedgeenable}" 
	    puts $datei "stloptions.resthedgeanglefac  ${stloptions.resthedgeanglefac}" 
	    puts $datei "stloptions.resthedgeangleenable  ${stloptions.resthedgeangleenable}" 
	    puts $datei "stloptions.resthsurfmeshcurvfac  ${stloptions.resthsurfmeshcurvfac}" 
	    puts $datei "stloptions.resthsurfmeshcurvenable  ${stloptions.resthsurfmeshcurvenable}" 
	    puts $datei "stloptions.recalchopt  ${stloptions.recalchopt}" 
	    
	    puts $datei "visoptions.subdivisions ${visoptions.subdivisions}"
	    puts $datei "visoptions.autoredraw ${visoptions.autoredraw}"
	    puts $datei "visoptions.autoredrawtime ${visoptions.autoredrawtime}"


	    # trafo options   
	    # if exist trafooptions then ...
	    if { [info exists trafooptions.solver] == 1 } {
		puts $datei "trafooptions.solver ${trafooptions.solver}" 
		puts $datei "trafooptions.levels ${trafooptions.levels}" 
		puts $datei "trafooptions.linits ${trafooptions.linits}" 
		puts $datei "trafooptions.nonlinits ${trafooptions.nonlinits}" 
		puts $datei "trafooptions.stabcurrent ${trafooptions.stabcurrent}" 
		puts $datei "trafooptions.checkcond ${trafooptions.checkcond}" 
		puts $datei "trafooptions.maxdirect ${trafooptions.maxdirect}" 
		puts $datei "trafooptions.secondorder ${trafooptions.secondorder}" 
		puts $datei "trafooptions.homogenizedcore ${trafooptions.homogenizedcore}" 
		puts $datei "trafooptions.ordercore ${trafooptions.ordercore}" 
		puts $datei "trafooptions.simplecurrents ${trafooptions.simplecurrents}" 
		puts $datei "trafooptions.assemblecomplexmatrix ${trafooptions.assemblecomplexmatrix}" 

		puts $datei "trafooptions.meshcasing  ${trafooptions.meshcasing}" 
		puts $datei "trafooptions.meshcore    ${trafooptions.meshcore}" 
		puts $datei "trafooptions.meshclumps  ${trafooptions.meshclumps}" 
		puts $datei "trafooptions.meshshields ${trafooptions.meshshields}" 
		puts $datei "trafooptions.meshcoils   ${trafooptions.meshcoils}" 
		puts $datei "trafooptions.bcmdirectory  ${trafooptions.bcmdirectory}" 
		puts $datei "trafooptions.lossdensityfile  ${trafooptions.lossdensityfile}" 
	    }

	    if { [info exists smalltrafomodell.tankheight] == 1 } {
		puts $datei "smalltrafomodell.tankheight ${smalltrafomodell.tankheight}"
		puts $datei "smalltrafomodell.tankwidth ${smalltrafomodell.tankwidth}"
		puts $datei "smalltrafomodell.tanklength ${smalltrafomodell.tanklength}"
		puts $datei "smalltrafomodell.corewidth ${smalltrafomodell.corewidth}"
		puts $datei "smalltrafomodell.windowheight ${smalltrafomodell.windowheight}"
		puts $datei "smalltrafomodell.limbdistance ${smalltrafomodell.limbdistance}"
		puts $datei "smalltrafomodell.xposcore ${smalltrafomodell.xposcore}"
		puts $datei "smalltrafomodell.yposcore ${smalltrafomodell.yposcore}"
		puts $datei "smalltrafomodell.zposcore ${smalltrafomodell.zposcore}"
		puts $datei "smalltrafomodell.leakagefluxguidethickness ${smalltrafomodell.leakagefluxguidethickness}"
		puts $datei "smalltrafomodell.leakagefluxguidewidth ${smalltrafomodell.leakagefluxguidewidth}"
		puts $datei "smalltrafomodell.leakagefluxguidezposition ${smalltrafomodell.leakagefluxguidezposition}"
		puts $datei "smalltrafomodell.limbcoil.1 ${smalltrafomodell.limbcoil.1}"
		puts $datei "smalltrafomodell.ricoil.1 ${smalltrafomodell.ricoil.1}"
		puts $datei "smalltrafomodell.rocoil.1 ${smalltrafomodell.rocoil.1}"
		puts $datei "smalltrafomodell.zposcoil.1 ${smalltrafomodell.zposcoil.1}"
		puts $datei "smalltrafomodell.heightcoil.1 ${smalltrafomodell.heightcoil.1}"
		puts $datei "smalltrafomodell.currentcoil.1 ${smalltrafomodell.currentcoil.1}"
		puts $datei "smalltrafomodell.nturnscoil.1 ${smalltrafomodell.nturnscoil.1}"
		puts $datei "smalltrafomodell.limbcoil.2 ${smalltrafomodell.limbcoil.2}"
		puts $datei "smalltrafomodell.ricoil.2 ${smalltrafomodell.ricoil.2}"
		puts $datei "smalltrafomodell.rocoil.2 ${smalltrafomodell.rocoil.2}"
		puts $datei "smalltrafomodell.zposcoil.2 ${smalltrafomodell.zposcoil.2}"
		puts $datei "smalltrafomodell.heightcoil.2 ${smalltrafomodell.heightcoil.2}"
		puts $datei "smalltrafomodell.currentcoil.2 ${smalltrafomodell.currentcoil.2}"
		puts $datei "smalltrafomodell.nturnscoil.2 ${smalltrafomodell.nturnscoil.2}"
		puts $datei "smalltrafomodell.limbcoil.3 ${smalltrafomodell.limbcoil.3}"
		puts $datei "smalltrafomodell.ricoil.3 ${smalltrafomodell.ricoil.3}"
		puts $datei "smalltrafomodell.rocoil.3 ${smalltrafomodell.rocoil.3}"
		puts $datei "smalltrafomodell.zposcoil.3 ${smalltrafomodell.zposcoil.3}"
		puts $datei "smalltrafomodell.heightcoil.3 ${smalltrafomodell.heightcoil.3}"
		puts $datei "smalltrafomodell.currentcoil.3 ${smalltrafomodell.currentcoil.3}"
		puts $datei "smalltrafomodell.nturnscoil.3 ${smalltrafomodell.nturnscoil.3}"
		puts $datei "smalltrafomodell.limbcoil.4 ${smalltrafomodell.limbcoil.4}"
		puts $datei "smalltrafomodell.ricoil.4 ${smalltrafomodell.ricoil.4}"
		puts $datei "smalltrafomodell.rocoil.4 ${smalltrafomodell.rocoil.4}"
		puts $datei "smalltrafomodell.zposcoil.4 ${smalltrafomodell.zposcoil.4}"
		puts $datei "smalltrafomodell.heightcoil.4 ${smalltrafomodell.heightcoil.4}"
		puts $datei "smalltrafomodell.currentcoil.4 ${smalltrafomodell.currentcoil.4}"
		puts $datei "smalltrafomodell.nturnscoil.4 ${smalltrafomodell.nturnscoil.4}"
		puts $datei "smalltrafomodell.limbcoil.5 ${smalltrafomodell.limbcoil.5}"
		puts $datei "smalltrafomodell.ricoil.5 ${smalltrafomodell.ricoil.5}"
		puts $datei "smalltrafomodell.rocoil.5 ${smalltrafomodell.rocoil.5}"
		puts $datei "smalltrafomodell.zposcoil.5 ${smalltrafomodell.zposcoil.5}"
		puts $datei "smalltrafomodell.heightcoil.5 ${smalltrafomodell.heightcoil.5}"
		puts $datei "smalltrafomodell.currentcoil.5 ${smalltrafomodell.currentcoil.5}"
		puts $datei "smalltrafomodell.nturnscoil.5 ${smalltrafomodell.nturnscoil.5}"
		puts $datei "smalltrafomodell.limbcoil.6 ${smalltrafomodell.limbcoil.6}"
		puts $datei "smalltrafomodell.ricoil.6 ${smalltrafomodell.ricoil.6}"
		puts $datei "smalltrafomodell.rocoil.6 ${smalltrafomodell.rocoil.6}"
		puts $datei "smalltrafomodell.zposcoil.6 ${smalltrafomodell.zposcoil.6}"
		puts $datei "smalltrafomodell.heightcoil.6 ${smalltrafomodell.heightcoil.6}"
		puts $datei "smalltrafomodell.currentcoil.6 ${smalltrafomodell.currentcoil.6}"
		puts $datei "smalltrafomodell.nturnscoil.6 ${smalltrafomodell.nturnscoil.6}"
		puts $datei "smalltrafomodell.limbtest.1 ${smalltrafomodell.limbtest.1}"
		puts $datei "smalltrafomodell.heighttest.1 ${smalltrafomodell.heighttest.1}"
		puts $datei "smalltrafomodell.widthtest.1 ${smalltrafomodell.widthtest.1}"
		puts $datei "smalltrafomodell.rtest.1 ${smalltrafomodell.rtest.1}"
		puts $datei "smalltrafomodell.zpostest.1 ${smalltrafomodell.zpostest.1}"
		puts $datei "smalltrafomodell.edgeradiustest.1 ${smalltrafomodell.edgeradiustest.1}"
		puts $datei "smalltrafomodell.finetest.1 ${smalltrafomodell.finetest.1}"
		puts $datei "smalltrafomodell.conductivetest.1 ${smalltrafomodell.conductivetest.1}"
		puts $datei "smalltrafomodell.limbtest.2 ${smalltrafomodell.limbtest.2}"
		puts $datei "smalltrafomodell.heighttest.2 ${smalltrafomodell.heighttest.2}"
		puts $datei "smalltrafomodell.widthtest.2 ${smalltrafomodell.widthtest.2}"
		puts $datei "smalltrafomodell.rtest.2 ${smalltrafomodell.rtest.2}"
		puts $datei "smalltrafomodell.zpostest.2 ${smalltrafomodell.zpostest.2}"
		puts $datei "smalltrafomodell.edgeradiustest.2 ${smalltrafomodell.edgeradiustest.2}"
		puts $datei "smalltrafomodell.finetest.2 ${smalltrafomodell.finetest.2}"
		puts $datei "smalltrafomodell.conductivetest.2 ${smalltrafomodell.conductivetest.2}"
		puts $datei "smalltrafomodell.limbtest.3 ${smalltrafomodell.limbtest.3}"
		puts $datei "smalltrafomodell.heighttest.3 ${smalltrafomodell.heighttest.3}"
		puts $datei "smalltrafomodell.widthtest.3 ${smalltrafomodell.widthtest.3}"
		puts $datei "smalltrafomodell.rtest.3 ${smalltrafomodell.rtest.3}"
		puts $datei "smalltrafomodell.zpostest.3 ${smalltrafomodell.zpostest.3}"
		puts $datei "smalltrafomodell.edgeradiustest.3 ${smalltrafomodell.edgeradiustest.3}"
		puts $datei "smalltrafomodell.finetest.3 ${smalltrafomodell.finetest.3}"
		puts $datei "smalltrafomodell.conductivetest.3 ${smalltrafomodell.conductivetest.3}"
		puts $datei "smalltrafomodell.limbtest.4 ${smalltrafomodell.limbtest.4}"
		puts $datei "smalltrafomodell.heighttest.4 ${smalltrafomodell.heighttest.4}"
		puts $datei "smalltrafomodell.widthtest.4 ${smalltrafomodell.widthtest.4}"
		puts $datei "smalltrafomodell.rtest.4 ${smalltrafomodell.rtest.4}"
		puts $datei "smalltrafomodell.zpostest.4 ${smalltrafomodell.zpostest.4}"
		puts $datei "smalltrafomodell.edgeradiustest.4 ${smalltrafomodell.edgeradiustest.4}"
		puts $datei "smalltrafomodell.finetest.4 ${smalltrafomodell.finetest.4}"
		puts $datei "smalltrafomodell.conductivetest.4 ${smalltrafomodell.conductivetest.4}"
		puts $datei "smalltrafomodell.nperitest ${smalltrafomodell.nperitest}"
		puts $datei "smalltrafomodell.filename ${smalltrafomodell.filename}"
		puts $datei "smalltrafomodell.murlfguide ${smalltrafomodell.murlfguide}"
		puts $datei "smalltrafomodell.murtestwire ${smalltrafomodell.murtestwire}"
		puts $datei "smalltrafomodell.murcore ${smalltrafomodell.murcore}"
		puts $datei "smalltrafomodell.kappalfguide ${smalltrafomodell.kappalfguide}"
		puts $datei "smalltrafomodell.kappatestwire ${smalltrafomodell.kappatestwire}"
		puts $datei "smalltrafomodell.kappacore ${smalltrafomodell.kappacore}"
	    }
	    
	    
	    close $datei
	}
    }
}




# the ini file is saved on demand :
proc saveinifile { } {
    global inifilename
    if {[catch { set datei [open $inifilename w] } result ]} {
	puts "cannot write file $inifilename"
    } {
	for { set i [.ngmenu.file.recent index last] } { $i >= 0 } { incr i -1 } {
	    puts $datei "recentfile \"[.ngmenu.file.recent entrycget $i -label]\""
	}
	close $datei
    }    
}


proc savemeshinifile { } {
    global meshinifilename 
    if {[catch { set datei [open $meshinifilename w] } result ]} {
	puts "cannot write file $meshinifilename"
    } {
	for { set i [.ngmenu.file.recentmesh index last] } { $i >= 1 } { incr i -1 } {
	    puts $datei "recentfile \"[.ngmenu.file.recentmesh entrycget $i -label]\""
	}
	close $datei
    }    
}



proc loadinifile { } { 
    global inifilename
    if { [file exists $inifilename] == 1 } {
	set datei [open $inifilename r]
	while { [gets $datei line] >= 0 } {
	    if {[lindex $line 0] == "recentfile"} {
		    set filename [lindex $line 1]
		    if { [file exists $filename] == 1 } {
		        AddRecentFile $filename
		    }	
	    }
	}
	close $datei
    }
}


proc loadmeshinifile { } {
    global meshinifilename
    if { [file exists $meshinifilename] == 1 } {
	set datei [open $meshinifilename r]
	while { [gets $datei line] >= 0 } {
	    if {[lindex $line 0] == "recentfile"} {
		set filename [lindex $line 1]
		if { [file exists $filename] == 1 } {
		    AddRecentMeshFile $filename
		}	
	    }
	}
	close $datei
    }
 }





set cmdindex {}
set hlpindex {}
set secindex {}
