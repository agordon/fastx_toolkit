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
#ifndef __FASTX_ARGS__
#define __FASTX_ARGS__

#ifdef __cplusplus
extern "C" {
#endif

//One day this would all be OO :-)

const char* get_input_filename();
const char* get_output_filename();
int verbose_flag();
int compress_output_flag();
int get_fastq_ascii_quality_offset();
FILE* get_report_file();

typedef int (*parse_argument_func)(int optind, int optc, char* optarg)  ;

int fastx_parse_cmdline( int argc, char* argv[],
			 const char* program_options,
			 parse_argument_func program_parse_arg ) ;


#ifdef __cplusplus
}
#endif

#endif

