AC_INIT([ngsolve],[4.9.5],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX

AC_PREFIX_DEFAULT(["/opt/netgen"])
AC_DISABLE_STATIC
AC_PROG_LIBTOOL


AC_ARG_WITH([netgen], 
        [  --with-netgen=dir       use Netgen installed in directory dir],
        [netgendir=$withval],
        [netgendir="/opt/netgen"]
        )

AC_SUBST([NETGEN_INCLUDES], ["-DNETGEN_ELTRANS -I$netgendir/include"])

AC_CONFIG_HEADERS(config.h) AC_CONFIG_FILES(Makefile include/Makefile
ngstd/Makefile basiclinalg/Makefile linalg/Makefile fem/Makefile
comp/Makefile solve/Makefile multigrid/Makefile parallel/Makefile 
programming_demos/Makefile )

AC_OUTPUT
