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
"usage: fastx_reverse_complement [-h] [-r] [-z] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

FASTX fastx;

char reverse_complement_base ( const char input ) 
{
	switch(input)
	{
	case 'N':
		return 'N';
	case 'n':
		return 'n';
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	case 'a':
		return 't';
	case 't':
		return 'a';
	case 'g':
		return 'c';
	case 'c':
		return 'g';
	default:
		errx(1,"Invalid nucleotide value (%c) in reverse_complement_base()", input ); 
	}

	return '0'; //should not get here - just to please the compiler
}

void reverse_complement_fastx(FASTX* pFASTX)
{
	int i,j ;
	int length = strlen(pFASTX->nucleotides);

	char temp_nuc;
	int  temp_qual;

	for (i=0;i<length;i++)
		pFASTX->nucleotides[i] = reverse_complement_base ( pFASTX->nucleotides[i] ) ;

	i = 0 ;
	j = length - 1 ;
	while ( i < j ) {
		//Swap the nucleotides
		temp_nuc = pFASTX->nucleotides[i] ;
		pFASTX->nucleotides[i] = pFASTX->nucleotides[j] ;
		pFASTX->nucleotides[j] = temp_nuc;

		//Swap the quality scores
		if (pFASTX->read_fastq) {
			temp_qual = pFASTX->quality[i];
			pFASTX->quality[i] = pFASTX->quality[j];
			pFASTX->quality[j] = temp_qual ;
		}
		
		//Advance to next position
		i++;
		j--;
	}
}


int main(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "", NULL);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {
		reverse_complement_fastx(&fastx);
		fastx_write_record(&fastx);
	}

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Printing Reverse-Complement Sequences.\n" );
		fprintf(get_report_file(), "Input: %zu reads.\n", num_input_reads(&fastx) ) ;
		fprintf(get_report_file(), "Output: %zu reads.\n", num_output_reads(&fastx) ) ;
	}
	return 0;
}
