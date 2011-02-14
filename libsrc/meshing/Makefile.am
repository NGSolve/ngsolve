AM_CPPFLAGS = $(MPI_INCLUDES) -I$(top_srcdir)/libsrc/include 


noinst_HEADERS = adfront2.hpp hpref_quad.hpp meshfunc.hpp ruler3.hpp  \
adfront3.hpp findip.hpp findip2.hpp hpref_segm.hpp meshing2.hpp	      \
specials.hpp bisect.hpp geomsearch.hpp hpref_tet.hpp meshing3.hpp     \
topology.hpp boundarylayer.hpp global.hpp hpref_trig.hpp meshing.hpp  \
validate.hpp classifyhpel.hpp hpref_hex.hpp improve2.hpp meshtool.hpp \
clusters.hpp hprefinement.hpp improve3.hpp meshtype.hpp		      \
 hpref_prism.hpp localh.hpp msghandler.hpp curvedelems.hpp	      \
 hpref_pyramid.hpp meshclass.hpp ruler2.hpp bcfunctions.hpp	      \
 basegeom.hpp 



METASOURCES = AUTO

lib_LTLIBRARIES = libmesh.la

libmesh_la_SOURCES = adfront2.cpp adfront3.cpp bisect.cpp boundarylayer.cpp \
	clusters.cpp curvedelems.cpp delaunay.cpp delaunay2d.cpp	    \
	geomsearch.cpp global.cpp hprefinement.cpp improve2.cpp		    \
	improve2gen.cpp improve3.cpp localh.cpp meshclass.cpp		    \
	meshfunc.cpp meshfunc2d.cpp meshing2.cpp meshing3.cpp		    \
	meshtool.cpp meshtype.cpp msghandler.cpp netrule2.cpp		    \
	netrule3.cpp parser2.cpp parser3.cpp prism2rls.cpp		    \
	pyramid2rls.cpp pyramidrls.cpp quadrls.cpp refine.cpp		    \
	ruler2.cpp ruler3.cpp secondorder.cpp smoothing2.5.cpp		    \
	smoothing2.cpp smoothing3.cpp specials.cpp tetrarls.cpp		    \
	topology.cpp triarls.cpp validate.cpp zrefine.cpp bcfunctions.cpp   \
	parallelmesh.cpp  paralleltop.cpp  paralleltop.hpp basegeom.cpp 
