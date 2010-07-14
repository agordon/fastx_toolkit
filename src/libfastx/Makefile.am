
noinst_LIBRARIES = libfastx.a

AM_CPPFLAGS = $(CC_WARNINGS)

INCLUDES = $(GTEXTUTILS_CFLAGS) -I../libfastx
LDADD = $(GTEXTUTILS_LIBS) ../libfastx/libfastx.a $(LT_LDFLAGS)

libfastx_a_SOURCES = chomp.c chomp.h \
		     fastx.c fastx.h \
		     fastx_args.c fastx_args.h \
		     sequence_alignment.h sequence_alignment.cpp \
		     sequence.h sequence.cpp \
		     file_type_detector.h file_type_detector.cpp \
		     sequence_list_loader.h sequence_list_loader.cpp \
		     sequence_list_container.h sequence_list_container.cpp \
		     fasta_file.h fasta_file.cpp \
		     fastq_file.h fastq_file.cpp \
		     fastx_file.h fastx_file.cpp \
		     tab_file.h tab_file.cpp \
		     commandline_parameters.h commandline_parameters.cpp
