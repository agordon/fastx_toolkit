/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
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
"usage: fastq_quality_trimmer [-h] [-v] [-t N] [-l N] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-t N]       = Quality threshold - nucleotides with lower \n" \
"                  quality will be trimmed (from the end of the sequence).\n" \
"   [-l N]       = Minimum length - sequences shorter than this (after trimming)\n" \
"                  will be discarded. Default = 0 = no minimum length. \n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTQ input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTQ output file. default is STDOUT.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"\n";

int min_quality_threshold=0;
int min_length=0;

FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
	case 'l':
		if (optarg==NULL) 
			errx(1, "[-l] parameter requires an argument value");
		min_length = strtoul(optarg,NULL,10);
		if (min_length<0)
			errx(1,"Invalid minimum length value (-l %s)", optarg);
		break;

	case 't':
		if (optarg==NULL) 
			errx(1, "[-t] parameter requires an argument value");
		min_quality_threshold = strtol(optarg,NULL,10);
		break;
	default:
		errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	int i ;

	fastx_parse_cmdline(argc, argv, "t:l:", parse_program_args);

	if ( min_quality_threshold == 0 )
		errx(1, "Missing minimum quality threshold value (-t)" ) ;

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTQ_ONLY, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		//Scan each sequence - backwards
		for ( i=(int)strlen(fastx.nucleotides)-1 ; i >=0 ; i-- ) {
			if ( fastx.quality[i] < min_quality_threshold ) 
				fastx.nucleotides[i] = 0 ;
			else
				break ;	
		}

		if ( i>=0 && i+1 >= min_length )
			fastx_write_record(&fastx);
	}
	
	//
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Minimum Quality Threshold: %d\n", min_quality_threshold);
		if ( min_length > 0 )
			fprintf(get_report_file(), "Minimum Length: %d\n", min_length);
		else
			fprintf(get_report_file(), "No minimum Length\n");


		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;

		size_t discarded = num_input_reads(&fastx) - num_output_reads(&fastx) ;
		fprintf(get_report_file(), "discarded %zu (%zu%%) too-short reads.\n", 
			discarded, (discarded*100)/( num_input_reads(&fastx) ) ) ;
	}	

	return 0;
}
