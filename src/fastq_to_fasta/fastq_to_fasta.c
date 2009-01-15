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
"usage: fastq_to_fasta [-h] [-r] [-n] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
"\n" \
"version " VERSION "\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-r]         = Rename sequence identifiers to numbers.\n" \
"   [-n]         = keep sequences with unknown (N) nucleotides.\n" \
"                  Default is to discard such sequences.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA output file. default is STDOUT.\n" \
"\n";

FASTX fastx;
int flag_rename_seqid = 0;
int flag_discard_N = 1 ;

int parse_program_args(int optind, int optc, char* optarg)
{
	switch(optc) {
	case 'n':
		flag_discard_N = 0 ;
		break;
	
	case 'r':
		flag_rename_seqid = 1;
		break;
	}
	return 1;
}


int main(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "rn", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTQ_ONLY, ALLOW_N, REQUIRE_UPPERCASE);

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_FASTA, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

#ifdef DEBUG
		fastx_debug_print_record($fastx);
#endif

		//See if the input sequence contained 'N' nucleotides
		if ( flag_discard_N  && (strchr(fastx.nucleotides,'N') != NULL)) 
				continue;

		if ( flag_rename_seqid ) 
			snprintf(fastx.name, sizeof(fastx.name), "%zu", num_output_reads(&fastx)+1) ;

		fastx_write_record(&fastx);
	}

	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;

		if ( flag_discard_N ) {
			size_t discarded = num_input_reads(&fastx) - num_output_reads(&fastx) ;
			fprintf(get_report_file(), "discarded %u (%u%%) low-quality reads.\n", 
				discarded, (discarded*100)/( num_input_reads(&fastx) ) ) ;
		}
	}	

	return 0;
}
