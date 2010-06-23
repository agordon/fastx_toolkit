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
#include <climits>
#include <error.h>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <istream>
#include <getopt.h>
#include <cstring>
#include <err.h>
#include <errno.h>

#include <string>
#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "libfastx/sequence.h"
#include "libfastx/fastx_file.h"
#include "libfastx/tab_file.h"

using namespace std;
using namespace std::tr1;

string input_filename;
string output_filename;
bool verbose = false;
bool tabular_input = false;
bool tabular_output = false;

enum LONG_OPTIONS {
	OPT_TABULAR_INPUT = CHAR_MAX+1,
	OPT_TABULAR_OUTPUT
};

struct option filter_options[] = {
	{"tabin",	0,	NULL,	OPT_TABULAR_INPUT},
	{"tabout",	0,	NULL,	OPT_TABULAR_OUTPUT},
	{NULL,0,0,0},
};

void show_help()
{
	cout << "fastx_reader_writer_test [-v] [-h] [--tabin] [--tabout] [-i INPUT] [-o OUTPUT] [INPUT-FILES]" << endl;
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;
	int option_index;

	while ( (c=getopt_long(argc,argv,"i:o:hv", filter_options, &option_index )) != -1 ) {
		switch(c)
		{
			/* Standard Options */
		case 'i':
			input_filename = optarg;
			break;

		case 'h':
			show_help();
			break;

		case 'v':
			verbose = true ;
			break ;

		case 'o':
			output_filename = optarg;
			break;

		case OPT_TABULAR_INPUT:
			tabular_input = true;
			break;

		case OPT_TABULAR_OUTPUT:
			tabular_output = true;
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
}

int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);

	parse_command_line(argc, argv);

	ISequenceReader *pReader=NULL;
	ISequenceWriter *pWriter=NULL;

	if (tabular_input) {
		pReader = new TabularFileReader(input_filename, 64);
	} else {
		pReader = create_fastx_reader(input_filename, 64);
	}
	if (tabular_output) {
		pWriter = pReader->create_tabular_writer(output_filename);
	} else {
		pWriter = pReader->create_fastx_writer(output_filename);
	}

	Sequence seq;
	while (pReader->read_next_sequence(seq)) {
		pWriter->write_sequence(seq);
	}

	return 0;
}
