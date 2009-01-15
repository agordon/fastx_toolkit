#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <err.h>

#include <config.h>

#include "fastx.h"
#include "fastx_args.h"

const char* usage=
"usage: fastq_qual_conv [-h] [-a] [-n] [-z] [-i INFILE] [-f OUTFILE]\n" \
"\n" \
"version " VERSION "\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-a]         = Output ASCII quality scores (default).\n" \
"   [-n]         = Output numeric quality scores.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA output file. default is STDOUT.\n" \
"\n";

FASTX fastx;
int flag_output_ascii = 1;

int parse_program_args(int optind, int optc, char* optarg)
{
	switch(optc) {
	case 'a': //this is the default, nothing to change
		break;
	
	case 'n':
		flag_output_ascii = 0 ;
		break;
	}
	return 1;
}


int main(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "an", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTQ_ONLY, ALLOW_N, REQUIRE_UPPERCASE);

	fastx_init_writer(&fastx, get_output_filename(), 
		flag_output_ascii ? OUTPUT_FASTQ_ASCII_QUAL : OUTPUT_FASTQ_NUMERIC_QUAL,
		compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {
#ifdef DEBUG
		fastx_debug_print_record($fastx);
#endif
		fastx_write_record(&fastx);
	}

	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;
	}
	return 0;
}
