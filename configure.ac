# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([FASTX Toolkit],
        [0.0.13],
        [Assaf Gordon gordon@cshl.edu],
        [fastx_toolkit])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL
AX_C_LONG_LONG
AX_CXX_HEADER_STDCXX_TR1
if test "$ax_cv_cxx_stdcxx_tr1" != yes; then
	AC_MSG_ERROR([Your version of gcc does not support the 'std::tr1' standard. Recommended gcc version is 4.1.2 or later. Please use a newer gcc version, or try to download the pre-compiled binaries from the fastx-toolkit website.])
fi

PKG_CHECK_MODULES([GTEXTUTILS],[gtextutils])

dnl --enable-wall
EXTRA_CHECKS="-Wall -Wextra -Wformat-nonliteral -Wformat-security -Wswitch-default -Wswitch-enum -Wunused-parameter -Wfloat-equal -Werror"
AC_ARG_ENABLE(wall,
[  --enable-wall          Enable many common GCC warnings (-Wall,-Wextra, -Werror etc., default enabled)],
[case "${enableval}" in
  yes) wall=true ;;
  no)  wall=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-wall) ;;
esac],[wall=true])
if test "$wall" = "true"
then
  CFLAGS="${CFLAGS} ${EXTRA_CHECKS}"
  CXXFLAGS="${CXXFLAGS} ${EXTRA_CHECKS}"
fi

dnl --enable-debug
AC_ARG_ENABLE(debug,
[  --enable-debug          Enable debug mode (default enabled)],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=true])
if test "$debug" = "true"
then
  CFLAGS="${CFLAGS} -DDEBUG -g -O1"
  CXXFLAGS="${CFLAGS} -DDEBUG -g -O1"
else
  CFLAGS="${CFLAGS} -O3"
  CXXFLAGS="${CFLAGS} -O3"
fi

dnl 'all-static' marco copied from subversion's configure.ac
dnl Check for --enable-all-static option
AC_ARG_ENABLE(all-static,
	      AS_HELP_STRING([--enable-all-static],
			     [Build completely static (standalone) binaries.]),
			     [
			      if test "$enableval" = "yes" ; then
				      LT_LDFLAGS="-all-static $LT_LDFLAGS"
			      elif test "$enableval" != "no" ; then
				      AC_MSG_ERROR([--enable-all-static doesn't accept argument])
				      fi
				      ])

AC_SUBST(LT_LDFLAGS)



AC_CONFIG_FILES([
   Makefile
   doc/Makefile
   m4/Makefile
   src/Makefile
   src/libfastx/Makefile
   src/fastx_clipper/Makefile
   src/fastq_to_fasta/Makefile
   src/fastx_quality_stats/Makefile
   src/fastq_quality_converter/Makefile
   src/fastx_trimmer/Makefile
   src/fastq_quality_filter/Makefile
   src/fastq_quality_trimmer/Makefile
   src/fastx_artifacts_filter/Makefile
   src/fastx_reverse_complement/Makefile
   src/fastx_collapser/Makefile
   src/fastx_uncollapser/Makefile
   src/seqalign_test/Makefile
   src/fasta_formatter/Makefile
   src/fasta_nucleotide_changer/Makefile
   src/fastx_renamer/Makefile
   src/fastq_masker/Makefile
   galaxy/Makefile
   galaxy/tools/Makefile
   galaxy/tools/fastx_toolkit/Makefile
   galaxy/tools/cshl_custom_fastx_toolkit/Makefile
   galaxy/test-data/Makefile
   galaxy/static/Makefile
   galaxy/static/fastx_icons/Makefile
   galaxy/tool-data/Makefile
   scripts/Makefile
   build_scripts/Makefile
])

AC_OUTPUT
