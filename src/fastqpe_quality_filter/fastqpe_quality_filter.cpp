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
#include <iomanip>

#include <string>
#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "libfastx/sequence.h"
#include "libfastx/fastx.h"
#include "libfastx/fastx_file.h"
#include "libfastx/tab_file.h"

using namespace std;
using namespace std::tr1;

string input_filename[2];
string output_filename[2];
bool verbose = false;
bool tabular_input = true;
bool tabular_output = true;
int min_quality=20 ;
int min_percent=100;
int ASCII_quality_offset = 64 ;

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
"fastqpe_quality_filter [-v] [-h]\n" \
"                           [--out1 OUTPUT1] [--out2 OUTPUT2]\n" \
"                           [-q N] [-p N] [-Q N]\n" \
"                           [--in1 INPUT1] [--in2 INPUT2]\n" \
"                           [INPUT1] [INPUT2]\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-q N]       = Minimum quality score to keep.\n" \
"   [-p N]       = Minimum percent of bases that must have [-q] quality.\n" \
"   [--in1 INPUT1]  = FASTQ input file (First of the paired-end files)\n" \
"   [--in2 INPUT2]  = FASTQ input file (Second of the paired-end files)\n" \
"   [--out1 FILE1]  = FASTQ output file (first of the paired-end files)\n" \
"   [--out2 FILE2]  = FASTQ output file (second of the paired-end files)\n" \
"                     If output files are not specified, output will be sent to STDOUT\n" \
"                     As tabulated FASTQ files (8 columns per pair-sequences).\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [--out] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"\n";
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;
	char *endptr;
	int option_index;

	while ( (c=getopt_long(argc,argv,"q:p:hv", filter_options, &option_index )) != -1 ) {
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

		case 'Q':
			ASCII_quality_offset = (int)strtol(optarg, &endptr, 10);
			if (endptr == optarg) {
				cerr << "Parameter error: invalid ASCII-quality-offset value (-q " << optarg << ")" << endl;
				exit(1);
			}
			break;
		case 'q':
			min_quality = (int)strtol(optarg, &endptr, 10);
			if (endptr == optarg) {
				cerr << "Parameter error: invalid minimum quality value (-q " << optarg << ")" << endl;
				exit(1);
			}
			break;

		case 'p':
			min_percent = (int)strtol(optarg, &endptr, 10);
			if (endptr == optarg) {
				cerr << "Parameter error: invalid minimum percent value (-p " << optarg << ")" << endl;
				exit(1);
			}
			if (min_percent<0 || min_percent>100) {
				cerr << "Parameter error: invalid minimum percent value (-p " << optarg << "), must be a number between 0 and 100." << endl;
				exit(1);
			}
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

int get_index_of_nth_element(int *array, int array_size, int n)
{
	int pos;

	//Find the first nono-empty index
	pos = 0 ;
	while ( pos < array_size && array[pos]==0 )
		pos++;

	#if 0
	fprintf(stderr,"n=%d\n", n);
	for (i=0; i< array_size; i++) {
		if (array[i] != 0)
			fprintf(stderr, "[%d]=%d  ", i + MIN_QUALITY_VALUE, array[i]) ;
	}
	fprintf(stderr,"\n");
	#endif
	
	if (pos == array_size)
		errx(1,"bug: got empty array at %s:%d", __FILE__, __LINE__);
	
	while (n > 0) {
		if (array[pos] > n)
			break;
		n -= array[pos];
		pos++;
		while (array[pos]==0 && pos < array_size)
			pos++;
	}
	return pos;
}

int get_percentile_quality(const Sequence& seq1, const Sequence& seq2, int percentile)
{
	size_t i;
	int count=0;
	int quality_values[QUALITY_VALUES_RANGE];

	memset(quality_values, 0, sizeof(quality_values));

	for (i=0; i< seq1.quality_cached_line.length(); ++i) {
		count++;
		quality_values[ (int)seq1.quality_cached_line[i] - ASCII_quality_offset - MIN_QUALITY_VALUE ] ++ ;
	}
	for (i=0; i< seq2.quality_cached_line.length(); ++i) {
		count++;
		quality_values[ (int)seq2.quality_cached_line[i] - ASCII_quality_offset - MIN_QUALITY_VALUE ] ++ ;
	}

	i = get_index_of_nth_element(quality_values, QUALITY_VALUES_RANGE, (count * (100-percentile) / 100));
	
	//printf(" n = %d, i = %d, i+MIN_QUAL_VALUE=%d\n", 
	//	(count*(100-percentile)/100), i, i+MIN_QUALITY_VALUE) ;

	return i + MIN_QUALITY_VALUE ;
}

int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);
	size_t seq_count_in=0;
	size_t seq_count_out=0;

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
		++seq_count_in;
		if (! (seq1.ASCII_quality_scores && seq1.ASCII_quality_scores)) {
			cerr << "Input error: Numeric quality scores are not supported by this program. Please convert you FASTQ files to ASCII quality scores" << endl;
			exit(1);
		}

		int value = get_percentile_quality(seq1,seq2, min_percent);
		if (value >= min_quality) {
			++seq_count_out;
			pWriter->write_sequence(seq1, seq2);
		}
	}
	if (verbose) {
		ostream &os(tabular_output?cerr:cout);

		os << "Input:  " << seq_count_in << " sequences." << endl;
		os << "Output: " << seq_count_out << " sequences." << endl;
		os << "Filtered-out:" << (seq_count_in-seq_count_out) <<
			"(" << fixed << setw(3) << ( ((double)((seq_count_in-seq_count_out)*100))/(seq_count_in*1.0)) << "%)" << endl;
	}

	return 0;
}
