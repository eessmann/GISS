#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(FES_distribution_soft, 1.0.9, tletelli@notos.cst.cnes.fr)
AC_CONFIG_SRCDIR([src/fes2004_alloc.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE(FES2004_distribution_soft, 1.0.9)
# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

#CFLAGS="$CFLAGS -g"



######################
#                    #
# Checks for include #
#                    #
######################
dnl --------RECHERCHE DE include NETCDF--------
AC_CHECK_HEADER(netcdf.h,,AC_MSG_ERROR([

Can not found netcdf.h
install it or use >> export CFLAGS="\$CFLAGS -I/netcdf_PATH/include"
]))



######################
#                    #
# Checks for library #
#                    #
######################

dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB C--------
dnl 
dnl ----------------------------------------

AC_CHECK_LIB(c,main)



dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB pthread--------
dnl 
dnl ----------------------------------------

#AC_CHECK_LIB(pthread,pthread_create ,  , [AC_MSG_ERROR([can not found the library libpthread the software can be compiled but generally it will not work !!!!!])])
AC_CHECK_LIB(pthread,pthread_create)



dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB MATH--------
dnl 
dnl ----------------------------------------

AC_CHECK_LIB(m,main)


dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB BLAS--------
dnl 
dnl ----------------------------------------

AC_CHECK_LIB(gslcblas,main, ,[AC_MSG_ERROR([

Can not found libgslcblas
install it or use :
1- export LDFLAGS="\$LDFLAGS -L/gsl_PATH/lib"
2- export LDFLAGS="\$LDFLAGS -L$(gsl-config --prefix)/lib"

])])

dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB GSL--------
dnl 
dnl ----------------------------------------

AC_CHECK_LIB(gsl,main, ,[AC_MSG_ERROR([

Can not found libgsl
install it or use :
1- export LDFLAGS="\$LDFLAGS -L/gsl_PATH/lib"
2- export LDFLAGS="\$LDFLAGS -L$(gsl-config --prefix)/lib"

])])

dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB NETCDF--------
dnl 
dnl ----------------------------------------

AC_CHECK_LIB([netcdf],[main], ,[AC_MSG_ERROR([

Can not found libnetcdf
install it or use >> export LDFLAGS="\$LDFLAGS -L/netcdf_PATH/lib"

])])



AC_PROG_RANLIB

AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor sqrt strdup])

dnl ----------------------------------------
dnl 
dnl --------RECHERCHE DE LIB fes--------
dnl 
dnl ----------------------------------------
dnl La lib FES fait planter la recherche des autres composants....
dnl il faut faire le test en dernier pour s'affranchir des erreurs
dnl mais il faut la faire pour l'include dans la compil 

AC_CHECK_LIB(fes,extraction, ,LIBS=" -lfes $LIBS")



AC_OUTPUT(Makefile src/Makefile)
