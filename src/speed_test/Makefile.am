# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


noinst_PROGRAMS = is_nuc_string \
		  cmdline_param_test \
		  fastxse_param_test

AM_CPPFLAGS = $(CC_WARNINGS)

INCLUDES = $(GTEXTUTILS_CFLAGS) -I$(top_srcdir)/src
LDADD = $(GTEXTUTILS_LIBS) ../libfastx/libfastx.a $(LT_LDFLAGS)

is_nuc_string_SOURCES = is_nuc_string.cpp
cmdline_param_test_SOURCES = cmdline_param_test.cpp
fastxse_param_test_SOURCES = fastxse_param_test.cpp
