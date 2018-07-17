#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libgsm], [1.3], [jaideep777@gmail.com])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_PROG_AR

AC_MSG_NOTICE([Setting up libtool])
LT_INIT

# Checks for programs.
: ${CXXFLAGS=""}	# remove the default -O2 -g flags
AC_PROG_CXX

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])

# Checks for libraries.
#  Compile application with netcdf
AC_LANG_PUSH([C++])
NETCDF_ROOT=''
NETCDF_INCLUDE=''
NETCDF_LIBS=''
ENABLE_NETCDF=no
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<yes|no|directory> (default=no)],[location of NetCDF-C++ library (lib and include subdirs)])],
            [AS_CASE(["$with_netcdf"],
                     [no],[AC_MSG_CHECKING([for NetCDF-C++ library])
                           AC_MSG_ERROR([NetCDF-C++ library is required.])],
                     [yes],[AC_CHECK_HEADERS([netcdfcpp.h])
#                            AC_SEARCH_LIBS([nc_open],
#                                           [netcdf],
#                                           [AC_DEFINE([HAVE_LIBNETCDF],[1],[Define to 1 for NetCDF support])
#                                            ENABLE_NETCDF=yes],
#                                           [AC_MSG_ERROR([Could not link to NetCDF library])])
                            AC_CHECK_HEADERS([netcdfcpp.h])
                            AX_CXX_CHECK_LIB(netcdf_c++, [NcFile::is_valid( void ) const])
                            NETCDF_LIBS=" -lnetcdf_c++"
                          	],
                     [*],[AS_IF([test -d "$with_netcdf"],
                                [NETCDF_ROOT=$with_netcdf
                                 LDFLAGS="-L$NETCDF_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$NETCDF_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([netcdfcpp.h])
                                 AX_CXX_CHECK_LIB(netcdf_c++, [NcFile::is_valid( void ) const])
                                 NETCDF_LIBS=" -L$NETCDF_ROOT/lib -lnetcdf_c++"
                                 NETCDF_INCLUDE=" -I$NETCDF_ROOT/include"
								],
                                [AC_MSG_NOTICE([$with_netcdf is not a directory! NetCDF-C++ suppressed])])
                         ]
                      )
            ],
            [AC_MSG_CHECKING([for NetCDF-C++ library])
             AC_MSG_ERROR([NetCDF-C++ library is required.])]
           )

AC_LANG_POP([C++])



AC_SUBST([ENABLE_NETCDF])
AC_SUBST([NETCDF_ROOT])
AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDF_LIBS])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([LDFLAGS])


AC_OUTPUT


