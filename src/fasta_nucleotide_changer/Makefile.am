# Copyright (C) 2008 Assaf Gordon <gordon@cshl.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


bin_PROGRAMS = fasta_nucleotide_changer

AM_CPPFLAGS = \
	$(CC_WARNINGS) \
	-I../libfastx

fasta_nucleotide_changer_SOURCES = fasta_nucleotide_changer.c

fasta_nucleotide_changer_LDADD = ../libfastx/libfastx.a $(LT_LDFLAGS)

