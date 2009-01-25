AC_INIT([netgen],[4.9.5],[],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT(["/opt/netgen"])


# Tcl/Tk configuration:

TEA_INIT([3.6])

TEA_PATH_TCLCONFIG
TEA_LOAD_TCLCONFIG 

TEA_PATH_TKCONFIG
TEA_LOAD_TKCONFIG

TEA_PUBLIC_TCL_HEADERS




AC_ARG_ENABLE([occ], 
        [  --enable-occ            compile with OpenCascade geometry kernel],
        [occon=true 
        occfl="-DOCCGEOMETRY"
        occlfl="-lTKIGES -lTKBRep -lTKSTEP -lTKSTL" ],
        [occon=false] )

AC_ARG_WITH([occ], 
        [  --with-occ=dir          use OpenCascade installed in directory dir],
        [occdir=$withval]
        [occincfl="-I$withval/inc"]
        [occlibfl="-L$withval/lib"],
        )


AC_ARG_WITH([togl], 
        [  --with-togl=dir         directory containing libTogl1.7],
        [togldir=$withval]
        [togllibfl="-L$withval$"]
        )
#        [togllibfl="-L$(TK_BIN_DIR)/Togl1.7"]
#        [togllibfl="-L/usr/local/lib/Togl1.7"]


AC_ARG_ENABLE([nglib],
        [  --enable-nglib          generate shared library nglib],
        [ nglibon=true ],
        [ nglibon=false ])


AC_ARG_ENABLE([parallel],
        [  --enable-parallel       enable mpi parallelization],
        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -DPARALLEL  -I/usr/share/metis/Lib -DMETIS")
        AC_SUBST([MPI_LIBS], "-lmetis -L/opt/mpich/ch-p4/lib -lmpich")             
        ]
        )

#        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -DPARALLEL  -I/home/joachim/download/metis-4.0/Lib -DMETIS")
#        AC_SUBST([MPI_LIBS], "-L/home/joachim/download/metis-4.0 -lmetis -L/opt/mpich/ch-p4/lib64 -lmpich")             


# AC_HEADER_STDC
AC_CONFIG_HEADERS(config.h)

AC_DISABLE_STATIC

AC_LANG([C++])
AC_PROG_CXX
AC_PROG_LIBTOOL
LT_INIT


# AM_CONDITIONAL([OCC], [test x$occon = xtrue])
AM_CONDITIONAL([NGLIB], [test x$nglibon = xtrue])

# AC_SUBST([OCCFLAGS], [$occfl])
AC_SUBST([OCCFLAGS], ["$occfl $occincfl"])
AC_SUBST([OCCLIBS], ["$occlfl $occlibfl"])
AC_SUBST([TOGLLIBDIR], ["$togllibfl"])

AC_CHECK_HEADER(pthread.h)
AC_CHECK_HEADER([togl.h])

AC_CHECK_HEADER([iostream],[OCCFLAGS="$OCCFLAGS -DHAVE_IOSTREAM"])
AC_CHECK_HEADER([iostream.h],[OCC_FLAGS="$OCCFLAGS -DHAVE_IOSTREAM_H"])
AC_CHECK_HEADER([limits],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS"])
AC_CHECK_HEADER([limits.h],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS_H"])

AC_CHECK_LIB(pthread, pthread_create)
# AC_SEARCH_LIBS(Tcl_Init, [tcl8.5 tcl8.4])


AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([floor])
AC_CHECK_HEADERS([limits.h])


AC_CONFIG_FILES(Makefile libsrc/Makefile libsrc/csg/Makefile 
libsrc/general/Makefile  libsrc/geom2d/Makefile libsrc/gprim/Makefile
libsrc/include/Makefile libsrc/interface/Makefile 
libsrc/linalg/Makefile libsrc/meshing/Makefile libsrc/occ/Makefile 
  libsrc/stlgeom/Makefile
  libsrc/visualization/Makefile ng/Makefile nglib/Makefile
  tutorials/Makefile doc/Makefile)

AC_OUTPUT
