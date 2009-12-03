# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


bin_PROGRAMS = fastq_quality_filter

AM_CPPFLAGS = \
	$(CC_WARNINGS) \
	-I../libfastx

fastq_quality_filter_SOURCES = fastq_quality_filter.c

fastq_quality_filter_LDADD = ../libfastx/libfastx.a $(LT_LDFLAGS)

