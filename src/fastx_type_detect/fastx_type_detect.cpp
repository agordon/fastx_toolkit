/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2010  A. Gordon (gordon@cshl.edu)

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
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <getopt.h>
#include <tr1/memory>
#include <cstring>
#include <err.h>
#include <errno.h>

#include <gtextutils/generic_input_stream.h>

using namespace std;

typedef std::tr1::shared_ptr<std::istream> shared_istream_ptr;

string input_filename;
bool allow_fasta = false;
bool allow_fastq = false;
bool allow_compressed = true;
bool check_nucleotides = true;
bool verbose = false;
bool allow_multiline = false;

struct option type_detection_options[] = {
	{"fasta",	0,	NULL,	'F'},
	{"fastq",	0,	NULL,	'Q'},
	{"fastx",	0,	NULL,	'X'},
	{"no-gzip",	0,	NULL,	'G'},
	{"no-nt-check",	0,	NULL,	'N'},
	{"allow-multiline",	0,	NULL,	'M'},
	{"help",	0,	NULL,	'h'},
	{"verbose",	0,	NULL,	'v'},
	{NULL,0,0,0},
};

void show_help()
{
	cout << "" << endl;
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;
	int option_index;
	bool default_behaviour = true ;

	while ( (c=getopt_long(argc,argv,"vi:FQXGNMh", type_detection_options, &option_index)) != -1 ) {
		switch(c)
		{
			/* Standard Options */
		case 'i':
			input_filename = optarg;
			break;

		case 'h':
			show_help();
			break;

		case 'F':
			allow_fasta = true;
			default_behaviour = false;
			break;

		case 'Q':
			allow_fastq = true;
			default_behaviour = false;
			break;

		case 'X':
			allow_fasta = allow_fastq = true;
			default_behaviour = false;
			break;

		case 'G':
			allow_compressed = false;
			default_behaviour = false;
			break;

		case 'v':
			verbose = true ;
			default_behaviour = false;
			break ;

		case 'M':
			allow_multiline = true ;
			default_behaviour = false;
			break;

		case 'N':
			check_nucleotides = false;
			default_behaviour = false;
			break;

		default:
			exit(1);
			break;
		}
	}
	//if no file name specified with "-i" and there's an extra argument - assume it is the file name
	if ( input_filename.empty() && optind < argc ) {
		input_filename = argv[optind];
	}

	//If no command line options specific, use the defaults: allow all, verbose
	if ( default_behaviour ) {
		verbose = true ;
	}
}

int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);

	parse_command_line(argc, argv);

	generic_input_stream a( input_filename ) ;

	cerr << "Sorry, this program is not finished yet." << endl;
	exit(1);

	return 0;
}
