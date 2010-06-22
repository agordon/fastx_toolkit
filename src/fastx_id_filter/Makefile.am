# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = fastx_id_filter

AM_CPPFLAGS = $(CC_WARNINGS)

INCLUDES = $(GTEXTUTILS_CFLAGS) -I../libfastx
LDADD = $(GTEXTUTILS_LIBS) ../libfastx/libfastx.a $(LT_LDFLAGS)

fastx_id_filter_SOURCES = fastx_id_filter.cpp \
				file_type_detector.cpp file_type_detector.h \
				sequence_list_loader.cpp sequence_list_loader.h \
				sequence_list_container.cpp sequence_list_container.h

