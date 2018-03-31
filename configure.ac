#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(arcs, 1.0.1, rwarren@bcgsc.ca)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([Arcs/Arcs.cpp])
AC_CONFIG_HEADER([config.h])

AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Options to configure.

# htslib (samtools parsing library)
AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib=PATH],
[specify prefix directory for htslib (samtools parsing library)]))
if test "$with_htslib" -a -d "$with_htslib"; then
	htslib_cppflags="-isystem$with_htslib/include"
	if test -d "$with_htslib/lib64"; then
		htslib_ldflags="-L$with_htslib/lib64"
	else
		htslib_ldflags="-L$with_htslib/lib"
	fi
fi

# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
			[specify directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
	boost_cppflags="-isystem$with_boost"
fi

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
boost_ver=1.58.0
boost_ver_dir=1_58_0
AC_SUBST(CPPFLAGS,"-I$srcdir $htslib_cppflags $boost_cppflags $CPPFLAGS -isystem$srcdir/$boost_ver_dir")
AC_SUBST(LDFLAGS, "$htslib_ldflags $LDFLAGS")

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers
AC_CHECK_HEADERS([boost/graph/graph_traits.hpp])
AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp])
AC_CHECK_HEADERS([boost/graph/undirected_graph.hpp])
AC_CHECK_HEADERS([boost/graph/graphviz.hpp])
#AC_CHECK_HEADERS([zlib.h])

# Check for htslib
AC_CHECK_HEADERS([htslib/sam.h])
if test x"$ac_cv_header_htslib_sam_h" != x"yes"; then
	AC_MSG_ERROR([could not find dependency htslib (try --with-htslib=<path>)])
fi

# Check for Boost.
if test $ac_cv_header_boost_graph_graphviz_hpp != yes; then
	AC_MSG_ERROR([Requires the Boost C++ libraries, which may
	be downloaded from here: http://www.boost.org/users/download/
	The following commands will download and install Boost:
	cd $srcdir
	wget http://downloads.sourceforge.net/project/boost/boost/$boost_ver/$boost_ver_dir.tar.bz2
	tar jxf $boost_ver_dir.tar.bz2
	cd -])
fi

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror -std=c++0x')

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		Arcs/Makefile
		Common/Makefile
		DataLayer/Makefile
		DataStructures/Makefile
		DistanceEst/Makefile
		Graph/Makefile
		Test/Makefile
])
AC_OUTPUT
