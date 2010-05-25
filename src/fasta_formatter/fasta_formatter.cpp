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
#include <cstring>
#include <err.h>
#include <vector>

#include <gtextutils/stream_wrapper.h>
#include <gtextutils/text_line_reader.h>
#include <gtextutils/string_tokenize.h>

#include "sequence_writers.h"

#include "config.h"

using namespace std;

string input_filename;
string output_filename;
bool flag_output_empty_sequences = false ;
bool flag_output_tabular = false ;
int  flag_requested_output_width = 0 ;
bool flag_input_tabular = false ;
int  tabular_input_seqid_column = 1 ;
int  tabular_input_nucl_column = 2 ;

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
"   [-T]	 = Input is a Tabular file, not a fasta file.\n" \
"   [-T n]	     n = identifier column (defualt=1)\n" \
"   [-T n,m]	     m = sequence column (default=n+1)\n" \
"\n" \
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

void parse_tabular_column_argument(const char* optarg)
{
	char *endptr;

	tabular_input_seqid_column = (int)strtoul(optarg, &endptr, 10);
	if ( endptr == optarg )
		errx(1,"Error: invalid argument to -T '%s'. Must be a numeric value.", optarg);
	if ( tabular_input_seqid_column == 0 )
		errx(1,"Error: invalid argument to -T '%s'. Must be a numeric value > 0.", optarg);

	//If there are two columns specifed
	if ( *endptr == ',' ) {
		char *second_arg = endptr+1;
		tabular_input_nucl_column = (int)strtoul(second_arg, &endptr, 10);
		if ( endptr == optarg )
			errx(1,"Error: invalid argument to -T '%s'. Expecting two numeric values separated by a comma (e.g. -T 2,6).", optarg);
		if ( tabular_input_nucl_column == 0 )
			errx(1,"Error: invalid second argument to -T '%s'. Must be a numeric value > 0.", optarg);
	} else {
		tabular_input_nucl_column = tabular_input_seqid_column + 1 ;
	}
}

void parse_command_line(int argc, char* argv[])
{
	int opt;

	while ( (opt = getopt(argc, argv, "i:o:hw:teT::") ) != -1 ) {

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

		case 'T':
			flag_input_tabular = true ;
			//Check for column specification from the user (in form of "n[,m]"
			if ( optarg != NULL && strlen(optarg)>0)
				parse_tabular_column_argument(optarg) ;
			break;

		default:
			exit(1);
		}
	}
}

void reformat_fasta_input( TextLineReader &reader, SequencesWriter &writer )
{
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
				writer.write ( sequence_id, sequence_bases ) ;

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
	writer.write ( sequence_id, sequence_bases ) ;
}

void reformat_tabular_input( TextLineReader &reader, SequencesWriter &writer )
{
	string sequence_id ;
	string sequence_bases ;
	while ( reader.next_line() ) {

		const string &line = reader.line_string();
		if ( line.length()==0 )
			continue;

		std::vector<string> fields;
		String_Tokenize ( line, std::back_inserter(fields), "\t" ) ;

		if (fields.size() >= (size_t)tabular_input_seqid_column)
			sequence_id = string(">") + fields[tabular_input_seqid_column-1];
		else
			errx(1,"Input error: line %zu is too short (expecting %d columns, for sequence-id field)\n",
					reader.line_number(), tabular_input_seqid_column ) ;

		if (fields.size() >= (size_t)tabular_input_nucl_column)
			sequence_bases = fields[tabular_input_nucl_column-1];
		else
			errx(1,"Input error: line %zu is too short (expecting %d columns, for nucleotides field)\n",
					reader.line_number(), tabular_input_nucl_column ) ;

		writer.write ( sequence_id, sequence_bases ) ;
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

	if ( flag_input_tabular ) {
		reformat_tabular_input ( reader, *pWriter ) ;
	} else {
		reformat_fasta_input ( reader, *pWriter ) ;
	}



	delete pWriter;
}

