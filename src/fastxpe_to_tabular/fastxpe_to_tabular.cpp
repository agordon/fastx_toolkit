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

string input_filename[2];
string output_filename[2];
bool verbose = false;
bool input_1_specified=false;
bool input_2_specified=false;

enum LONG_OPTIONS {
	OPT_INPUT_END1 = CHAR_MAX+1,
	OPT_INPUT_END2
};

struct option filter_options[] = {
	{"in1",		1,	NULL,	OPT_INPUT_END1},
	{"in2",		1,	NULL,	OPT_INPUT_END2},
	{NULL,0,0,0},
};

void show_help()
{
	cout <<
"fastxpe_to_tabular [-v] [-h]\n" \
"                           [--in1 INPUT1] [--in2 INPUT2]\n" \
"			or\n" \
"                           [INPUT1] [INPUT2]\n" \
"\n";
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;
	int option_index;

	while ( (c=getopt_long(argc,argv,"hv", filter_options, &option_index )) != -1 ) {
		switch(c)
		{
			/* Standard Options */
		case OPT_INPUT_END1:
			input_filename[0] = optarg;
			input_1_specified=true;
			break;

		case OPT_INPUT_END2:
			input_filename[1] = optarg;
			input_2_specified=true;
			break;

		case 'h':
			show_help();
			break;

		case 'v':
			verbose = true ;
			break ;

		default:
			exit(1);
			break;
		}
	}

	if ( input_1_specified != input_2_specified ) {
			cerr << "Parameters error: Only one input file specified. When reading from files, two input files are required (one for each end of the paired-end data)."  << endl;
			exit(1);
	}

	//if no file name specified with "--inX" and there are two extra arguments
        //  - assume they are two input file names
	if ( !input_1_specified ) {
		if (optind+1 < argc) {
			input_filename[0] = argv[optind];
			input_filename[1] = argv[optind+1];
		} else
		if (optind < argc) {
			cerr << "Parameters error: a single filename was given ("
				<< argv[optind] << "). Please provide two file names (one for each end of the paired-end data)." << endl;
			exit(1);
		} else {
			cerr << "Parameters error: No input files specified." << endl;
			exit(1);
		}
	}
}

int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);

	parse_command_line(argc, argv);

	ISequenceReaderPE *pReader=NULL;
	ISequenceWriterPE *pWriter=NULL;

	pReader = create_fastx_pe_reader(input_filename[0], input_filename[1], 64);
	pWriter = pReader->create_tabular_writer("");

	Sequence seq1;
	Sequence seq2;
	while (pReader->read_next_sequence(seq1, seq2)) {
		pWriter->write_sequence(seq1, seq2);
	}

	return 0;
}
