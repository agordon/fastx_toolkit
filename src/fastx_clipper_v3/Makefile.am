# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


bin_PROGRAMS = fastx_clipper_v3

AM_CPPFLAGS = \
	$(CC_WARNINGS) \
	-I../libfastx

INCLUDES = $(GTEXTUTILS_CFLAGS) -I../libfastx
LDADD = $(GTEXTUTILS_LIBS) ../libfastx/libfastx.a $(LT_LDFLAGS)

fastx_clipper_v3_SOURCES = fastx_clipper_v3.cpp \
			   kmer_coding.h kmer_coding.cpp \
			   adapter_hash.h adapter_hash.cpp \
			   sequence_writer.h \
			   tile_detection.h tile_detection.cpp \
			   fastx_writer.h fastx_writer.cpp

#			   sequence_alignment.h sequence_alignment.cpp \
#
