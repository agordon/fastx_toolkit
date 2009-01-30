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
        [0.0.4],
        [Assaf Gordon gordon@cshl.edu],
        [fastx_toolkit])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])

# 23dec08, Gordon
# Only added those things because 'autoheader' was complaining...
AC_DEFINE([CXX_HAS_BUGGY_FOR_LOOPS], [], [Description])
AC_DEFINE([CXX_HAS_NO_BOOL], [], [Description])
AC_DEFINE([NDEBUG], [], [Description])
AC_DEFINE([YOUR_OS], [], [Description])

#Is this the right way to add compiler flags?
CFLAGS="-Wall -g"
CXXFLAGS="-Wall -g"


LF_CONFIGURE_CC
LF_CONFIGURE_CXX
LF_HOST_TYPE
LF_SET_WARNINGS
AC_PROG_RANLIB

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
   src/fastx_artifacts_filter/Makefile
   src/fastx_reverse_complement/Makefile
   src/fastx_collapser/Makefile
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
