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
"usage: fastq_masker [-h] [-v] [-q N] [-r C] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-q N]       = Quality threshold - nucleotides with lower quality will be masked\n" \
"                  Default is 10.\n" \
"   [-r C]       = Replace low-quality nucleotides with character C. Default is 'N'\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTQ input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTQ output file. default is STDOUT.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"\n";

int min_quality_threshold=10;
char mask_character='N';

FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
	case 'q':
		if (optarg==NULL)
			errx(1, "[-q] parameter requires an argument value");
		min_quality_threshold = atoi(optarg);
		if (min_quality_threshold<-40)
			errx(1,"Invalid minimum length value (-q %s)", optarg);
		break;

	case 'r':
		if (optarg==NULL)
			errx(1, "[-r] parameter requires an argument value");
		if (strlen(optarg)!=1)
			errx(1, "[-r] paramter requires a single character as value");
		mask_character = optarg[0];
		break;

	default:
		errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	int i ;
	size_t masked_reads_count=0;
	size_t masked_nucleotides_count=0;

	fastx_parse_cmdline(argc, argv, "q:r:", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(),
		FASTQ_ONLY, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		int masked = 0;

		//Scan each sequence - backwards
		for ( i=0; i<(int)strlen(fastx.nucleotides); ++i ) {
			if ( fastx.quality[i] < min_quality_threshold ) {
				fastx.nucleotides[i] = mask_character ;
				masked = 1;
				++masked_nucleotides_count;
			}
		}
		if (masked)
			masked_reads_count += get_reads_count(&fastx);

		fastx_write_record(&fastx);
	}
	//
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Minimum Quality Threshold: %d\n", min_quality_threshold);
		fprintf(get_report_file(), "Low-quality nucleotides replaced with '%c'\n", mask_character);

		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;

		fprintf(get_report_file(), "Masked reads: %zu\n", masked_reads_count ) ;
		fprintf(get_report_file(), "Masked nucleotides: %zu\n", masked_nucleotides_count ) ;
	}

	return 0;
}
