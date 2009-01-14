AC_INIT([netgen], [4.9.2], [], [])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


AC_ARG_ENABLE([occ], 
        [  --enable-occ            Compile with OpenCascade geometry kernel],
        [occon=true occfl='-DOCCGEOMETRY -I/usr/local/inc'
        occlfl='-lTKIGES -lTKBRep -lTKSTEP -lTKSTL' ],
        [occon=false] )


AC_CONFIG_HEADERS(config.h)

AC_LANG_CPLUSPLUS
AC_PROG_CXX
AM_PROG_LIBTOOL

AC_SYS_POSIX_TERMIOS

AM_CONDITIONAL([OCC], [test x$occon = xtrue])
AC_SUBST([OCCFLAGS], [$occfl])
AC_SUBST([OCCLIBS], [$occlfl])

AC_CHECK_HEADER([iostream],[OCCFLAGS="$OCCFLAGS -DHAVE_IOSTREAM"])
AC_CHECK_HEADER([iostream.h],[OCC_FLAGS="$OCCFLAGS -DHAVE_IOSTREAM_H"])
AC_CHECK_HEADER([limits],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS"])
AC_CHECK_HEADER([limits.h],[OCCFLAGS="$OCCFLAGS -DHAVE_LIMITS_H"])



AC_CONFIG_FILES(Makefile libsrc/Makefile libsrc/csg/Makefile libsrc/general/Makefile \
libsrc/geom2d/Makefile libsrc/gprim/Makefile libsrc/include/Makefile libsrc/interface/Makefile \
libsrc/linalg/Makefile libsrc/meshing/Makefile libsrc/occ/Makefile libsrc/opti/Makefile \
libsrc/parallel/Makefile libsrc/stlgeom/Makefile libsrc/visualization/Makefile ng/Makefile \
tutorials/Makefile)

AC_OUTPUT



