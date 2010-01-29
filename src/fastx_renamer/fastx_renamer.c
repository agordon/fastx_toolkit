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
"usage: fastx_renamer [-n TYPE] [-h] [-z] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-n TYPE]    = rename type:\n" \
"                  SEQ - use the nucleotides sequence as the name.\n" \
"                  COUNT - use simply counter as the name.\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

enum RENAME_TYPE {
	SEQ,
	COUNT
}; 
enum RENAME_TYPE rename_type;
unsigned int counter = 1 ;
FASTX fastx;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
	case 'n':
		if (optarg==NULL) 
			errx(1, "[-n] parameter requires an argument value");
		if (strncmp(optarg,"SEQ",3)==0) {
			rename_type = SEQ ;
		}
		else
		if (strncmp(optarg,"COUNT",5)==0) {
			rename_type = COUNT ;
			counter = 1 ;
		}
		else 
			errx(1,"Uknown rename type [-n]: '%s'", optarg);
		break;

	default:
		errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;

	}
	return 1;
}


int main(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "n:", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE, get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		switch(rename_type)
		{
		case SEQ:
			strncpy(fastx.name, fastx.nucleotides, sizeof(fastx.name));
			strncpy(fastx.name2, fastx.nucleotides, sizeof(fastx.name2));
			break;
		case COUNT:
			snprintf(fastx.name, sizeof(fastx.name),"%u",counter);
			strncpy(fastx.name2, fastx.name, sizeof(fastx.name2));
			counter++;
			break;
		default:
			errx(1,"Internal error: rename_type = %d", (int)rename_type);
		}

		fastx_write_record(&fastx);
	}

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Renamed: %zu reads.\n", num_input_reads(&fastx) ) ;
	}
	return 0;
}
