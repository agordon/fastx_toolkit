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
"usage: fastx_trimmer [-h] [-f N] [-l N] [-z] [-v] [-i INFILE] [-o OUTFILE]\n" \
"\n" \
"version " VERSION "\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-f N]       = First base to keep. Default is 1 (=first base).\n" \
"   [-l N]       = Last base to keep. Default is entire read.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

#define DO_NOT_TRIM_LAST_BASE (0)

//Default adapter (Solexa's MODBAN) used at the Hannon Lab
int keep_first_base=1;
int keep_last_base=DO_NOT_TRIM_LAST_BASE;

FASTX fastx;

int parse_program_args(int optind, int optc, char* optarg)
{
	switch(optc) {
	case 'f':
		if (optarg==NULL) 
			errx(1, "[-f] parameter requires an argument value");
		keep_first_base = strtoul(optarg,NULL,10);
		if (keep_first_base<=0 || keep_first_base>=MAX_SEQ_LINE_LENGTH) 
			errx(1,"Invalid number bases to keep (-f %s)", optarg);
		break;

	case 'l':
		if (optarg==NULL) 
			errx(1, "[-l] parameter requires an argument value");
		keep_last_base = strtoul(optarg,NULL,10);
		if (keep_last_base<=0 ||  keep_last_base>=MAX_SEQ_LINE_LENGTH) 
			errx(1,"Invalid number bases to keep (-l %s)", optarg);
		break;

	}
	return 1;
}


int main(int argc, char* argv[])
{
	int i;
	
	fastx_parse_cmdline(argc, argv, "l:f:", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE);

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		if (keep_last_base != DO_NOT_TRIM_LAST_BASE) {
			fastx.nucleotides[keep_last_base] = 0 ;
		}

		if (keep_first_base != 1) {
			for (i=0; i < strlen(fastx.nucleotides)-keep_first_base+1 ; i++) {
				fastx.nucleotides[i] = fastx.nucleotides[i+keep_first_base-1];
				fastx.quality[i] = fastx.quality[i+keep_first_base-1];
			}
			fastx.nucleotides[i] = 0 ;
		}

		//none of the above condition matched, so print this sequence.
		fastx_write_record(&fastx);
	}

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Trimming: base %d to %d\n", keep_first_base, keep_last_base ) ;
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;
	}
	return 0;
}
