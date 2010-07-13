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
bool tabular_input = true;
bool tabular_output = true;

enum LONG_OPTIONS {
	OPT_INPUT_END1 = CHAR_MAX+1,
	OPT_INPUT_END2,
	OPT_OUTPUT_END1,
	OPT_OUTPUT_END2
};

struct option filter_options[] = {
	{"in1",		1,	NULL,	OPT_INPUT_END1},
	{"in2",		1,	NULL,	OPT_INPUT_END2},
	{"out1",	1,	NULL,	OPT_OUTPUT_END1},
	{"out2",	1,	NULL,	OPT_OUTPUT_END2},
	{NULL,0,0,0},
};

void show_help()
{
	cout <<
"fastxpe_reader_writer_test [-v] [-h]\n" \
"                           [--in1 INPUT1] [--in2 INPUT2]\n" \
"                           [--out1 OUTPUT1] [--out2 OUTPUT2]\n" \
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
			tabular_input=false;
			break;

		case OPT_INPUT_END2:
			input_filename[1] = optarg;
			tabular_input=false;
			break;

		case OPT_OUTPUT_END1:
			output_filename[0] = optarg;
			tabular_output=false;
			break;

		case OPT_OUTPUT_END2:
			output_filename[1] = optarg;
			tabular_output=false;
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

	//if no file name specified with "--inX" and there are two extra arguments
        //  - assume they are two input file names
	if ( tabular_input ) {
		if (optind+1 < argc) {
			input_filename[0] = argv[optind];
			input_filename[1] = argv[optind+1];
			tabular_input = false;
		} else
		if (optind < argc) {
			cerr << "Parameters error: a single filename was given ("
				<< argv[optind] << "). Please provide two file names (one for each end of the paired-end data), or no file names at all (to read from STDIN)." << endl;
			exit(1);
		}
	}

	//validate possible bad combination of parameters
	if ( !tabular_input ) {
		//must have two input file names
		if (input_filename[0].empty() != input_filename[1].empty()) {
			cerr << "Parameters error: Only one input file specified. When reading from files, two input files are required (one for each end of the paired-end data)."  << endl;
			exit(1);
		}
	}
	if ( !tabular_output) {
		//must have two input file names
		if (output_filename[0].empty() != output_filename[1].empty()) {
			cerr << "Parameters error: Only one output file specified. When writing to files, two output files are required (one for each end of the paired-end data)."  << endl;
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

	if (tabular_input) {
		pReader = new PE_TabularFileReader("", 64);
	} else {
		pReader = create_fastx_pe_reader(input_filename[0], input_filename[1], 64);
	}
	if (tabular_output) {
		pWriter = pReader->create_tabular_writer("");
	} else {
		pWriter = pReader->create_fastx_writer(output_filename[0], output_filename[1]);
	}

	Sequence seq1;
	Sequence seq2;
	while (pReader->read_next_sequence(seq1, seq2)) {
		pWriter->write_sequence(seq1, seq2);
	}

	return 0;
}
