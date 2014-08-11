AC_INIT([ngsolve],[5.3-dev],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign nostdinc])

AC_CONFIG_MACRO_DIR([m4])

AC_LANG([C++])

AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

# from libtool
# LDFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"   


AC_PREFIX_DEFAULT(["/opt/netgen"])
AC_DISABLE_STATIC
AM_PROG_AR
AC_PROG_LIBTOOL



CXXFLAGS="$CXXFLAGS -DUSE_TIMEOFDAY -DLAPACK"
# ACX_LAPACK
                

AC_ARG_WITH([netgen], 
        [  --with-netgen=dir       use Netgen installed in directory dir],
        [netgendir=$withval],
        [netgendir="/opt/netgen"]
        )



guion=true
AC_ARG_ENABLE([gui],
        [AC_HELP_STRING([--disable-gui],[don't build ngsolve with GUI])],
        [if test "$enableval" = yes; then guion=true; else guion=false; fi])

AM_CONDITIONAL([NGSGUI], [test x$guion = xtrue])

if test a$guion = atrue ; then
# Tcl configuration:
TEA_INIT([3.9])
TEA_PATH_TCLCONFIG
TEA_LOAD_TCLCONFIG 
TEA_PUBLIC_TCL_HEADERS
fi

		

AC_ARG_WITH([lapack],
        [AC_HELP_STRING([--with-lapack=lapacklib][linker flags for linking dense linear algebra library LAPACK, default is '-llapack -lblas'])],          
        [AC_SUBST([LAPACK_LIBS], [$withval])],
        [AC_SUBST([LAPACK_LIBS], "-llapack -lblas")],
        )

AC_ARG_ENABLE([pardiso],
        [AC_HELP_STRING([--enable-pardiso],[enable sparse direct solver pardiso])], 
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO -DUSE_PARDISO400"])]
        [AC_SUBST([PARDISO_LIBS], ["-lpardiso500-GNU481-X86-64 -lgfortran"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )
#  -llapack

AC_ARG_ENABLE([pardiso4],
        [AC_HELP_STRING([--enable-pardiso4],[enable sparse direct solver pardiso])],
        [AC_SUBST([PARDISO_INCLUDES], ["$PARDISO_INCLUDES -DUSE_PARDISO400"])]
        )

AC_ARG_WITH([pardiso],
        [AC_HELP_STRING([--with-pardiso=libpardiso.so],[enable sparse direct solver pardiso])],
        [AC_SUBST([PARDISO_INCLUDES], ["$PARDISO_INCLUDES -DUSE_PARDISO"])]
        [AC_SUBST([PARDISO_LIBS], ["$withval -llapack -lblas -lgomp -gfortran"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )

# AC_ARG_WITH([superlu],
#        [AC_HELP_STRING([--with-superlu=-I/usr/include/superlu],[enable sparse LU factorization library])],
#        [AC_SUBST([SUPERLU_INCLUDES], ["-DUSE_SUPERLU $withval"])]
#        [AC_SUBST([SUPERLU_LIBS], ["-lsuperlu $withval"])]
#        )


AC_ARG_ENABLE([MKLpardiso],
        [  --enable-MKLpardiso        enable sparse direct solver pardiso from MKL],          
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO -DUSE_MKL"])]
        )

AC_ARG_ENABLE([mumps],
        [AC_HELP_STRING([--enable-mumps],[enable sparse direct solver MUMPS])],
        [AC_SUBST([MUMPS_INCLUDES], ["-DUSE_MUMPS"])]
        [AC_SUBST([MUMPS_LIBS], ["-L/opt/mumps/lib -ldmumps -lzmumps -lmumps_common -lpord -L/home/joachim/Downloads/scotch_5.1.12_esmumps/lib -lptesmumps -lptscotch -lptscotcherr  /usr/local/lib/libscalapack-2.0.2.a  /usr/local/lib/blacs_MPI-LINUX-0.a  /usr/local/lib/blacsCinit_MPI-LINUX-0.a  /usr/local/lib/blacs_MPI-LINUX-0-fake.a -lparmetis  -lmetis -lgfortran -lmpi_f77"])])
#  -lmpi_mpifh        
#  -lmpi_f77
#        [AC_SUBST([MUMPS_LIBS], ["-L/usr/local/lib -ldmumps -lzmumps -lmumps_common -lpord -L/home/joachim/Downloads/scotch_5.1.12_esmumps/lib   -lscalapack -lmpi_f90 -lgfortran"])] 
# -lgfortran
#        [AC_SUBST([MUMPS_LIBS], ["-L/opt/mumps/lib -ldmumps -lzmumps -lmumps_common -lpord -L/home/joachim/Downloads/scotch_5.1.12_esmumps/lib -lptesmumps -lptscotch -lptscotcherr  /usr/local/lib/libscalapack-2.0.2.a -lparmetis  -lmetis -lgfortran -lmpi_f77"])]



# AC_ARG_WITH([pardiso],
#       [  --with-pardiso=libpardiso.so  enable sparse direct solver pardiso],          
#       [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
#       [AC_SUBST([PARDISO_LIBS], [$withval])]
#       )

#	[AC_SUBST([LAPACK_LIB], ["-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_core -lmkl_def -lmkl_sequential -lmkl_gf_lp64"])]
#        [AC_SUBST([LAPACK_LIB], ["-llapack"])]
#	[AC_SUBST([LAPACK_LIB], ["-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_sequential"])]


#AC_SUBST([NETGEN_INCLUDES], ["-I$netgendir/include"])
# AC_SUBST(NETGEN_INCLUDES, "-I$netgendir/include")
# AC_MSG_NOTICE([Netgen includes... $NETGEN_INCLUDES])

CXXFLAGS="$CXXFLAGS -I. -I$netgendir/include"

AC_ARG_ENABLE([parallel],
        [AC_HELP_STRING([--enable-parallel],[enable mpi parallelization])],
        [CXXFLAGS="$CXXFLAGS -DPARALLEL"]
        )       

AC_ARG_ENABLE([vt],
        [AC_HELP_STRING([--enable-vt],[enable vampirtrace])],
        [CXXFLAGS="$CXXFLAGS -vt:inst manual -DVTRACE"]
        )       

AC_ARG_ENABLE([cuda],
	[AC_HELP_STRING([--enable-cuda],[enable CUDA GPU support])],
	[CXXFLAGS+=" -DCUDA -I/usr/local/cuda-6.5/include"]
	[LDFLAGS+=" -L/usr/local/cuda-6.5/lib64 -lcudart -lcublas -lcusparse"]
	)


echo "creating compiler wrapper ngscxx"
echo "$CXX -I. $GOLD_INCLUDE $CXXFLAGS -L$netgendir/lib -fpic $LAPACK_INCLUDES $TCL_INCLUDES $MUMPS_INCLUDES $PARDISO_INCLUDES \$*  $LAPACK_LIBS -Wl,--rpath=$netgendir/lib" > ngscxx
chmod a+x ngscxx 

define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))dnl
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])

define([builddate], esyscmd([sh -c "date|tr -d '\n'"]))dnl
AC_DEFINE(BUILDDATE, "builddate", [build date])

AC_SUBST(GOLD_INCLUDE, $GOLD_INCLUDE)

AC_CONFIG_HEADERS(config.h) 
AC_CONFIG_FILES(Makefile include/Makefile
ngstd/Makefile basiclinalg/Makefile linalg/Makefile fem/Makefile
comp/Makefile solve/Makefile multigrid/Makefile parallel/Makefile 
pde_tutorial/Makefile pde_tutorial/pml/Makefile windows/Makefile)

AC_OUTPUT
