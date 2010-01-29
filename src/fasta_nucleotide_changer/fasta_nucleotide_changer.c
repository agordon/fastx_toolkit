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
"usage: fasta_nucleotide_changer [-h] [-z] [-v] [-i INFILE] [-o OUTFILE] [-r] [-d]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-v]         = Verbose mode. Prints a short summary.\n" \
"                  with [-o], summary is printed to STDOUT.\n" \
"                  Otherwise, summary is printed to STDERR.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"   [-r]         = DNA-to-RNA mode - change T's into U's.\n" \
"   [-d]         = RNA-to-DNA mode - change U's into T's.\n" \
"\n";

int flag_dna_mode = 0;
int flag_rna_mode = 0;

FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char __attribute__((unused))* optarg)
{
	switch(optc) {
	case 'd':
		flag_dna_mode = 1 ;
		break;

	case 'r':
		flag_rna_mode = 1 ;
		break;

	default:
		errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;
	}
	return 1;
}


int main(int argc, char* argv[])
{
	size_t i;
	char nuc_from ;
	char nuc_to;
	size_t changes_count=0 ;
	
	fastx_parse_cmdline(argc, argv, "rd", parse_program_args);

	if ( !flag_dna_mode && !flag_rna_mode )
		errx(1,"Please specify either RNA mode (-r) or DNA mode (-d)" );

	if ( flag_dna_mode && flag_rna_mode )
		errx(1,"RNA mode (-r) and DNA mode (-d) can not be used together." );

	if ( flag_dna_mode ) {
		nuc_from = 'U';
		nuc_to   = 'T';
	}
	if ( flag_rna_mode ) {
		nuc_from = 'T';
		nuc_to   = 'U';
	}

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N | ALLOW_U, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_FASTA, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		for (i=0; i < strlen(fastx.nucleotides) ; i++)  {
			if ( fastx.nucleotides[i] == nuc_to ) {
				errx(1,"Error: found '%c' nucleotide on line %lld. (input should not contain '%c' nucleotides in %s mode)",
					nuc_to, fastx.input_line_number, nuc_to,(flag_dna_mode)?"RNA-to-DNA":"DNA-to-RNA" );
			}
			if ( fastx.nucleotides[i] == nuc_from ) {
				fastx.nucleotides[i] = nuc_to ;
				changes_count++;
			}
		}
		
		fastx_write_record(&fastx);
	}

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Mode: %s\n", (flag_dna_mode)?"RNA-to-DNA":"DNA-to-RNA" ) ;
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;
		fprintf(get_report_file(), "Nucleotides changed: %zu\n", changes_count ) ;
	}
	return 0;
}
