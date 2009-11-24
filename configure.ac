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
        [0.0.10],
        [Assaf Gordon gordon@cshl.edu],
        [fastx_toolkit])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

# 23dec08, Gordon
# Only added those things because 'autoheader' was complaining...
AC_DEFINE([CXX_HAS_BUGGY_FOR_LOOPS], [], [Description])
AC_DEFINE([CXX_HAS_NO_BOOL], [], [Description])
AC_DEFINE([NDEBUG], [], [Description])
AC_DEFINE([YOUR_OS], [], [Description])

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
   galaxy/Makefile
   galaxy/tools/Makefile
   galaxy/tools/fastx_toolkit/Makefile
   galaxy/tools/fastx_toolkit_with_gzip_and_output_label/Makefile
   galaxy/test-data/Makefile
   galaxy/static/Makefile
   galaxy/static/fastx_icons/Makefile
   galaxy/tool-data/Makefile
   scripts/Makefile
])

AC_OUTPUT
