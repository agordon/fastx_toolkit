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

using namespace std;
using namespace std::tr1;

string input_filename;
string output_filename;
bool verbose = false;

void show_help()
{
	cout << "fastx_reader_writer_test [-v] [-h] [-i INPUT] [-o OUTPUT] [INPUT-FILES]" << endl;
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;

	while ( (c=getopt(argc,argv,"i:o:hv")) != -1 ) {
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

	ISequenceReader *pReader = create_fastx_reader(input_filename, 64);
	ISequenceWriter *pWriter = pReader->create_writer(output_filename);

	Sequence seq;
	while (pReader->read_next_sequence(seq)) {
		pWriter->write_sequence(seq);
	}

	return 0;
}

