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

#define MAX_ADAPTER_LEN 100

const char* usage=
"usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-q N]       = Minimum quality score to keep.\n" \
"   [-p N]       = Minimum percent of bases that must have [-q] quality.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"\n";

#define DO_NOT_TRIM_LAST_BASE (0)

int min_quality=0;
int min_percent=0;

FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
	case 'q':
		if (optarg==NULL) 
			errx(1, "[-q] parameter requires an argument value");
		min_quality = strtoul(optarg,NULL,10);
		break;

	case 'p':
		if (optarg==NULL) 
			errx(1, "[-l] parameter requires an argument value");
		min_percent = strtoul(optarg,NULL,10);
		if (min_percent<=0 ||  min_percent>100) 
			errx(1,"Invalid percent value (-p %s)", optarg);
		break;
	default:
		errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;
	}
	return 1;
}

int get_index_of_nth_element(int *array, int array_size, int n)
{
	int pos;

	//Find the first nono-empty index
	pos = 0 ;
	while ( pos < array_size && array[pos]==0 )
		pos++;

	#if 0
	fprintf(stderr,"n=%d\n", n);
	for (i=0; i< array_size; i++) {
		if (array[i] != 0)
			fprintf(stderr, "[%d]=%d  ", i + MIN_QUALITY_VALUE, array[i]) ;
	}
	fprintf(stderr,"\n");
	#endif
	
	if (pos == array_size)
		errx(1,"bug: got empty array at %s:%d", __FILE__, __LINE__);
	
	while (n > 0) {
		if (array[pos] > n)
			break;
		n -= array[pos];
		pos++;
		while (array[pos]==0 && pos < array_size)
			pos++;
	}
	return pos;
}

int get_percentile_quality(const FASTX *fastx, int percentile)
{
	size_t i;
	int count=0;
	int quality_values[QUALITY_VALUES_RANGE];

	memset(quality_values, 0, sizeof(quality_values));

	for (i=0; i< strlen(fastx->nucleotides); i++) {
		count++;
		quality_values[ fastx->quality[i] - MIN_QUALITY_VALUE ] ++ ;
	}

	i = get_index_of_nth_element(quality_values, QUALITY_VALUES_RANGE, (count * (100-percentile) / 100));
	
	//printf(" n = %d, i = %d, i+MIN_QUAL_VALUE=%d\n", 
	//	(count*(100-percentile)/100), i, i+MIN_QUALITY_VALUE) ;

	return i + MIN_QUALITY_VALUE ;
}

int main(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "q:p:", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTQ_ONLY, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {
		#if 0
		fprintf(stderr, "%s\n", fastx.nucleotides ) ;
		for (i=0; i<strlen(fastx.nucleotides); i++) {
			fprintf(stderr,"%d ", fastx.quality[i]);
		}
		fprintf(stderr,"\n");
		#endif

		int value = get_percentile_quality(&fastx, min_percent);

		//fprintf(stderr, "value = %d\n\n", value ) ;
			

		if (value >= min_quality) {
			fastx_write_record(&fastx);
		} else {
	//		fprintf(stderr, "%s\n", fastx.nucleotides ) ;
	//		fprintf(stderr, "value = %d\n", value ) ;
		}
	}
	
	//
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Quality cut-off: %d\n", min_quality);
		fprintf(get_report_file(), "Minimum percentage: %d\n", min_percent);

		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;

		size_t discarded = num_input_reads(&fastx) - num_output_reads(&fastx) ;
		fprintf(get_report_file(), "discarded %zu (%zu%%) low-quality reads.\n", 
			discarded, (discarded*100)/( num_input_reads(&fastx) ) ) ;
	}	

	return 0;
}
