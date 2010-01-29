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
"usage: fastx_artifacts_filter [-h] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-v]         = Verbose - report number of processed reads.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"\n";

#define DO_NOT_TRIM_LAST_BASE (0)

FASTX fastx;

int parse_commandline(int argc, char* argv[])
{
	return fastx_parse_cmdline(argc, argv, "", NULL);
}

int artifact_sequence(const FASTX *fastx)
{
	int n_count=0;
	int a_count=0;
	int c_count=0;
	int t_count=0;
	int g_count=0;
	int total_count=0;

	int max_allowed_different_bases = 3 ;

	int i=0;

	while (1) {
		if (fastx->nucleotides[i]==0)
			break;

		total_count++;
		switch(fastx->nucleotides[i])
		{
		case 'A':
			a_count++;
			break;
		case 'C':
			c_count++;
			break;
		case 'G':
			g_count++;
			break;
		case 'T':
			t_count++;
			break;
		case 'N':
			n_count++;
			break;
		default:
			errx(1, __FILE__":%d: invalid nucleotide value (%c) at position %d",
				__LINE__, fastx->nucleotides[i], i ) ;
		}
		i++;
	}

	//Rules for artifacts
	
	if ( a_count>=(total_count-max_allowed_different_bases) 
	     ||
	     c_count>=(total_count-max_allowed_different_bases)
	     ||
	     g_count>=(total_count-max_allowed_different_bases)
	     ||
	     t_count>=(total_count-max_allowed_different_bases)
	     )
	     return 1;
	 

	 return 0;
}

int main(int argc, char* argv[])
{
	parse_commandline(argc, argv);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), 
		OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {
		
		if ( artifact_sequence(&fastx)  ) {
		} else {
			fastx_write_record(&fastx);
		}
	}
	
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;

		size_t discarded = num_input_reads(&fastx) - num_output_reads(&fastx) ;
		fprintf(get_report_file(), "discarded %zu (%zu%%) artifact reads.\n", 
			discarded, (discarded*100)/( num_input_reads(&fastx) ) ) ;
	}	

	return 0;
}
