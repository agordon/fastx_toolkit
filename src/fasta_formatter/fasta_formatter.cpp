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
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

#include <getopt.h>
#include <err.h>

#include <gtextutils/stream_wrapper.h>
#include <gtextutils/text_line_reader.h>

#include "sequence_writers.h"

#include "config.h"

using namespace std;

string input_filename;
string output_filename;
bool flag_output_empty_sequences = false ;
bool flag_output_tabular = false ;
int  flag_requested_output_width = 0 ;

const char* usage_string=
"usage: fasta_formatter [-h] [-i INFILE] [-o OUTFILE] [-w N] [-t] [-e]\n" \
"Part of " PACKAGE_STRING " by gordon@cshl.edu\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"   [-w N]       = max. sequence line width for output FASTA file.\n" \
"                  When ZERO (the default), sequence lines will NOT be wrapped -\n" \
"                  all nucleotides of each sequences will appear on a single \n" \
"                  line (good for scripting).\n" \
"   [-t]         = Output tabulated format (instead of FASTA format).\n" \
"                  Sequence-Identifiers will be on first column,\n" \
"                  Nucleotides will appear on second column (as single line).\n" \
"   [-e]         = Output empty sequences (default is to discard them).\n" \
"                  Empty sequences are ones who have only a sequence identifier,\n" \
"                  but not actual nucleotides.\n" \
"\n" \
"Input Example:\n" \
"   >MY-ID\n" \
"   AAAAAGGGGG\n" \
"   CCCCCTTTTT\n" \
"   AGCTN\n" \
"\n" \
"Output example with unlimited line width [-w 0]:\n" \
"   >MY-ID\n" \
"   AAAAAGGGGGCCCCCTTTTTAGCTN\n" \
"\n" \
"Output example with max. line width=7 [-w 7]:\n" \
"   >MY-ID\n" \
"   AAAAAGG\n" \
"   GGGTTTT\n" \
"   TCCCCCA\n" \
"   GCTN\n" \
"\n" \
"Output example with tabular output [-t]:\n" \
"   MY-ID	AAAAAGGGGGCCCCCTTTTAGCTN\n" \
"\n" \
"example of empty sequence:\n" \
"(will be discarded unless [-e] is used)\n" \
"  >REGULAR-SEQUENCE-1\n" \
"  AAAGGGTTTCCC\n" \
"  >EMPTY-SEQUENCE\n" \
"  >REGULAR-SEQUENCE-2\n" \
"  AAGTAGTAGTAGTAGT\n" \
"  GTATTTTATAT\n" \
"\n" \
"\n";

void usage()
{
	printf("%s",usage_string);
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int opt;

	while ( (opt = getopt(argc, argv, "i:o:hw:te") ) != -1 ) {
		
		//Parse the default options
		switch(opt) {
		case 'h':
			usage();
		
		case 'i':
			input_filename = optarg;
			break;

		case 'o':
			output_filename = optarg;
			break;

		case 'w':
			flag_requested_output_width = atoi(optarg);
			if ( flag_requested_output_width < 0 )
				errx(1,"Invalid value (%s) for requested width [-w]", optarg);
			break ;

		case 't':
			flag_output_tabular = true ;
			break;

		case 'e':
			flag_output_empty_sequences = true;
			break;
			
		default:
			exit(1);
		}
	}
}


int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);

	parse_command_line(argc, argv);

	InputStreamWrapper input ( input_filename ) ;
	OutputStreamWrapper output ( output_filename );
	TextLineReader reader ( input.stream() ) ;

	/*
	 * Use the writer according to the user's request
	 */
	SequencesWriter * pWriter = NULL ;

	if ( flag_output_tabular ) {
		pWriter = new TabulatedFastaWriter ( output.stream() ) ;
	} else {
		if ( flag_requested_output_width == 0 )
			pWriter = new SingleLineFastaWriter ( output.stream() ) ;
		else 
			pWriter = new MultiLineFastaWriter ( output.stream(), flag_requested_output_width ) ;
	}
	if (!flag_output_empty_sequences) {
		EmptySequencesFilter *filter = new EmptySequencesFilter ( pWriter ) ;
		pWriter = filter ;
	}


	/*
	 * FASTA read/process/write loop
	 */
	int max_length = 0 ;
	string sequence_id ;
	string sequence_bases ;
	bool first_line = true ;
	while ( reader.next_line() ) {

		const string &line = reader.line_string();
		
		if ( line.length()==0 )
			continue;

		if ( line[0] == '>' ) {
			//Got new sequence identifier - print previous sequence
			if (first_line)
				first_line = false;
			else
				pWriter->write ( sequence_id, sequence_bases ) ;
			
			// Start new sequence 
			sequence_id = line ;
			sequence_bases.clear();
			sequence_bases.resize ( max_length * 2 ) ;
		} else {
			//Got sequence nucleotides
			sequence_bases += line ;
		}
	}

	//Write the last sequence
	pWriter->write ( sequence_id, sequence_bases ) ;

	delete pWriter;
}

