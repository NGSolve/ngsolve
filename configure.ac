AC_INIT([netgen],[6.1-dev],[],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT(["/opt/netgen"])





# AC_HEADER_STDC
AC_CONFIG_HEADERS(config.h)
AC_DISABLE_STATIC

AC_LANG([C++])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
# LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"

AM_PROG_AR
AC_PROG_LIBTOOL
LT_INIT

# OpenCASCADE configuration

occon=false
AC_ARG_ENABLE([occ], 
        [AC_HELP_STRING([--enable-occ],[compile with OpenCascade geometry kernel])],
        [occon=true])

AC_ARG_WITH([occ], 
        [AC_HELP_STRING([--with-occ=dir],[use OpenCascade installed in directory dir])],
        [occdir=$withval] 
	[occon=true],
	[occdir=/opt/OpenCASCADE]
        )



if test a$occon = atrue ; then

	AC_SUBST([OCCFLAGS], ["-DOCCGEOMETRY -I$occdir/inc -I/usr/include/oce -I/usr/include/opencascade"])
	AC_SUBST([OCCLIBS], ["-L$occdir/lib -lTKernel -lTKGeomBase -lTKMath -lTKG2d -lTKG3d -lTKXSBase -lTKOffset -lTKFillet -lTKShHealing -lTKMesh -lTKMeshVS -lTKTopAlgo -lTKGeomAlgo -lTKBool -lTKPrim -lTKBO -lTKIGES -lTKBRep -lTKSTEPBase -lTKSTEP -lTKSTL -lTKSTEPAttr -lTKSTEP209 -lTKXDESTEP -lTKXDEIGES -lTKXCAF -lTKLCAF -lFWOSPlugin"])

#  -lTKDCAF

	if test a$build_cpu = ax86_64 ; then
		AC_SUBST([OCCFLAGS],["$OCCFLAGS -D_OCC64"])
	fi


	AC_CHECK_HEADER([iostream],[OCCFLAGS="$OCCFLAGS -DHAVE_IOSTREAM"])
	AC_CHECK_HEADER([iostream.h],[OCC_FLAGS="$OCCFLAGS -DHAVE_IOSTREAM_H"])
	AC_CHECK_HEADER([limits],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS"])
	AC_CHECK_HEADER([limits.h],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS_H"])
	AC_CHECK_HEADER([iomanip],[OCCFLAGS="$OCCFLAGS -DHAVE_IOMANIP"])
	AC_CHECK_HEADER([iomanip.h],[OCCFLAGS="$OCCFLAGS -DHAVE_IOMANIP_H"])


	echo "OCCFLAGS = $OCCFLAGS"
	echo "OCCLIBS = $OCCLIBS"
fi







AC_ARG_WITH([togl], 
        [AC_HELP_STRING([--with-togl=dir],[directory containing libTogl1.7])],
        [togldir=$withval]
        [togllibfl="-L$withval$"]
        )
#        [togllibfl="-L$(TK_BIN_DIR)/Togl1.7"]
#        [togllibfl="-L/usr/local/lib/Togl1.7"]


nglibon=true
AC_ARG_ENABLE([nglib],
        [AC_HELP_STRING([--enable-nglib],[generate shared library nglib])],
        [if test "$enableval" = yes; then nglibon=true; else nglibon=false; fi])

ngguion=true
AC_ARG_ENABLE([gui],
        [AC_HELP_STRING([--disable-gui],[don't build netgen with GUI])],
        [if test "$enableval" = yes; then ngguion=true; else ngguion=false; fi])

metisdir=/usr/local
AC_ARG_WITH([metis], 
        [AC_HELP_STRING([--with-metis=dir],[path to metis 5.x])],
        [metisdir=$withval],
        [metisdir=/usr/local]
	)

AC_ARG_ENABLE([parallel],
        [AC_HELP_STRING([--enable-parallel],[enable mpi parallelization])],
	[AC_SUBST([MPI_INCLUDES], "-I$metisdir/include -DMETIS")]
        [CXXFLAGS="$CXXFLAGS -DPARALLEL"]
	[AC_SUBST([MPI_LIBS], "-L$metisdir/lib -lmetis")]
        )       
#          -DVTRACE
#          -lvt-hyb

#        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -I/usr/include/mpich2 -I/usr/include/metis -DPARALLEL  -I/usr/share/metis/Lib -DMETIS -DVTRACE -I/usr/local/include/vampirtrace")
#        AC_SUBST([MPI_LIBS], "-lmetis -L/opt/mpich/ch-p4/lib -lmpich")

# -lvt -lvt-mpi -lvt-mpi-unify")
#        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -DPARALLEL -I/usr/include/metis -DMETIS")
#        AC_SUBST([MPI_LIBS], " -L/opt/mpich/ch-p4/lib64 -lmpich -lmetis")
#        [AC_SUBST([MPI_INCLUDES], "-DPARALLEL -DNO_PARALLEL_THREADS -I/opt/mpich/include")
#        AC_SUBST([MPI_LIBS], "")
# -lmetis
#        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -DPARALLEL  -I/home/joachim/download/metis-4.0/Lib -DMETIS")
#        AC_SUBST([MPI_LIBS], "-L/home/joachim/download/metis-4.0 -lmetis -L/opt/mpich/ch-p4/lib64 -lmpich")             


jpeglibon=false
AC_ARG_ENABLE([jpeglib],
        [AS_HELP_STRING([--enable-jpeglib],[enable snapshots using library libjpeg])],
        [AC_SUBST([JPEGLIB_INCLUDES], "-DJPEGLIB")
        AC_CHECK_LIB(jpeg, jpeg_start_compress,  
                     AC_SUBST([JPEGLIB_LIBS], "-ljpeg") , 
                     AC_MSG_ERROR([no usable library libjpeg found]))
        [if test "$enableval" = yes; then jpeglibon=true; else jpeglibon=false; fi]
        ]       
        )               

ffmpegon=false
AC_ARG_ENABLE([ffmpeg],
        [AS_HELP_STRING([--enable-ffmpeg],[enable video recording with FFmpeg, uses libavcodec])],
        [AC_SUBST([FFMPEG_INCLUDES], "-DFFMPEG -D__STDC_CONSTANT_MACROS")
        AC_SUBST([FFMPEG_LIBS], "-lavutil -lavformat -lavcodec -lavutil -lswscale -lz -lbz2")
        [if test "$enableval" = yes; then ffmpegon=true; else ffmpegon=false; fi]
        ]       
        )               


mklon=false
AC_ARG_ENABLE([mkl],
        [AS_HELP_STRING([--enable-mkl],[link Intel's mkl library, necessary for ngsolve and mkl 11.x])],
        [if test "$enableval" = yes; then mklon=true; else mklon=false; fi]
        )

AM_PATH_PYTHON([3.2])
AX_PYTHON_DEVEL([>= '3.2'])
AX_BOOST_BASE([1.54])
AX_BOOST_PYTHON

echo "pythondir = $pythondir"
# AC_SUBST([pythondir], \${prefix}/python/netgen)
AC_SUBST([pythondir], $pythondir/netgen)
# echo "pythondir changed to  = $pythondir"


pythonon=false
AC_ARG_ENABLE([python],
        [AS_HELP_STRING([--enable-python],[link boost python, necessary for ngsolve])],
        [if test "$enableval" = yes; then pythonon=true; else pythonon=false; fi]
	[CXXFLAGS+=" $PYTHON_CPPFLAGS -DNG_PYTHON"]
	[LDFLAGS+=" $PYTHON_LDFLAGS -l$BOOST_PYTHON_LIB"]
        )



# only for GUI version 

if test a$ngguion = atrue ; then
# Tcl/Tk configuration:
TEA_INIT([3.9])
TEA_PATH_TCLCONFIG
TEA_LOAD_TCLCONFIG 
TEA_PATH_TKCONFIG
TEA_LOAD_TKCONFIG
TEA_PUBLIC_TCL_HEADERS
TEA_PUBLIC_TK_HEADERS

AC_SUBST([TOGLLIBDIR], ["$togllibfl"])
# AC_CHECK_HEADER([togl.h])
AC_CHECK_HEADER([GL/gl.h])



#--------------------------------------------------------------------
# __CHANGE__
# Choose OpenGL platform (borrowed from Togl1.7)
#--------------------------------------------------------------------

case "${TEA_WINDOWINGSYSTEM}" in
	aqua)
		AC_SUBST(TOGL_WINDOWINGSYSTEM,TOGL_AGL)
		TEA_ADD_LIBS([-framework AGL -framework OpenGL -framework ApplicationServices])
		# libGLU is implicit in OpenGL framework
		LIBGLU=
		;;
	x11)
		AC_SUBST(TOGL_WINDOWINGSYSTEM,TOGL_X11)
		TEA_ADD_LIBS([-lGL -lXmu -lX11])
		LIBGLU=-lGLU
		;;
	win32)
		AC_SUBST(TOGL_WINDOWINGSYSTEM,TOGL_WGL)
		TEA_ADD_LIBS([opengl32.lib user32.lib gdi32.lib])
		if test "$GCC" = "yes" ; then
			LIBGLU=-lglu32
		else
			LIBGLU=glu32.lib
		fi
		;;
	*)
		AC_MSG_ERROR([Unsupported windowing system: ${TEA_WINDOWINGSYSTEM}])
		;;
esac
AC_SUBST(LIBGLU)

fi



AM_CONDITIONAL([NGLIB], [test x$nglibon = xtrue])
AM_CONDITIONAL([NGGUI], [test x$ngguion = xtrue])
AM_CONDITIONAL([NGMKL], [test x$mklon = xtrue])
AM_CONDITIONAL([NGPYTHON], [test x$pythonon = xtrue])

# AC_CHECK_HEADER(pthread.h)

# AC_CHECK_FUNCS([pow])
# AC_CHECK_FUNCS([floor])
# AC_CHECK_FUNCS([matherr])
# AC_CHECK_HEADERS([limits.h])

# AC_CHECK_LIB(pthread, pthread_create)
AC_CHECK_LIB(dl, dlsym)


# LDFLAGS+=" -L/usr/lib/Togl1.7"
# AC_SEARCH_LIBS(Togl_Init,[Togl Togl1.7],
#         MYTOGL_LIB=$LIBS
#        LIBS="",
#        ,
#        [-lGL -lXmu -ltk8.5 -ltcl8.5])
# echo "libs = $LIBS"
# echo "TOGL_LIB = $MYTOGL_LIB"

# AC_CONFIG_SUBDIRS([ng/Togl-1.7])

AC_CONFIG_FILES(Makefile libsrc/Makefile libsrc/csg/Makefile
libsrc/general/Makefile libsrc/geom2d/Makefile libsrc/gprim/Makefile
libsrc/include/Makefile libsrc/interface/Makefile
libsrc/linalg/Makefile libsrc/meshing/Makefile libsrc/occ/Makefile
libsrc/stlgeom/Makefile libsrc/visualization/Makefile ng/Makefile
nglib/Makefile tutorials/Makefile python/Makefile
py_tutorials/Makefile doc/Makefile windows/Makefile )

AC_OUTPUT








AC_MSG_RESULT([
------------------------------------------------------------------------
  $PACKAGE $VERSION:  Automatic configuration OK.

  Enabled functionality:

    OCC: ............... $occon
    JPEGlib: ........... $jpeglibon
    FFMPEG: ............ $ffmpegon
    NGLIB: ............. $nglibon
    GUI: ............... $ngguion


  Building:

    Type 'make' to compile $PACKAGE.

    Type 'make install' to install $PACKAGE.

  Example programs will be built but not installed.
------------------------------------------------------------------------
])
