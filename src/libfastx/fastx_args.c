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
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <getopt.h>

#include "fastx_args.h"

/*
 * Each program should specify its own usage string
 */
extern char* usage;


/*
 * globals.. yuck
 *
 * some day this will be a stand alone class
 */
const char* input_filename = "-";
const char* output_filename = "-";
int verbose = 0;
int compress_output = 0 ;
int fastq_ascii_quality_offset = 64 ;
FILE* report_file;

int get_fastq_ascii_quality_offset()
{
	return fastq_ascii_quality_offset;
}

const char* get_input_filename()
{
	return input_filename;
}

const char* get_output_filename()
{
	return output_filename;
}

int verbose_flag()
{
	return verbose;
}

int compress_output_flag()
{
	return compress_output ;
}

FILE* get_report_file()
{
	return report_file;
}

int fastx_parse_cmdline( int argc, char* argv[],
			 const char* program_options,
			 parse_argument_func program_parse_args ) 
{
	int opt;

	char combined_options_string[100];

	strcpy(combined_options_string, "Q:zhvi:o:");
	strcat(combined_options_string, program_options);
	
	report_file = stderr ; //since the default output is STDOUT, the report goes by default to STDERR

	while ( (opt = getopt(argc, argv, combined_options_string) ) != -1 ) {
		
		// Parse the program's custom options
		if ( strchr(program_options, opt) != NULL ) {
			if (!program_parse_args(optind, opt, optarg))
				return 0;
			continue;
		}

		//Parse the default options
		switch(opt) {
		case 'h':
			printf("%s", usage);
			exit(1);
		
		case 'v':
			verbose = 1 ;
			break ;

		case 'z':
			compress_output = 1 ;
			break ;


		case 'i':
			if (optarg==NULL)
				errx(1,"[-i] option requires FILENAME argument");
			input_filename = optarg;
			break;

		case 'o':
			if (optarg==NULL)
				errx(1,"[-o] option requires FILENAME argument");
			output_filename = optarg;
			
			//The user specified a specific output file, so the report can go to STDOUT
			report_file = stdout;
			break;
			
		case 'Q':
			if (optarg==NULL)
				errx(1,"[-Q] option requires VALUE argument");
			fastq_ascii_quality_offset = atoi(optarg);
			break;

		default:
			printf("use '-h' for usage information.\n");
			exit(1);
			break;

		}
	}

	return 1;
}

