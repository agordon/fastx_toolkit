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
#include <getopt.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <ios>
#include <iterator>
#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <fstream>

#include <gtextutils/stream_wrapper.h>
#include <gtextutils/text_line_reader.h>
#include <gtextutils/string_tokenize.h>
#include <gtextutils/container_join.h>

#include "config.h"

#include "fastx.h"
#include "fastx_args.h"

using namespace std;

const char* usage=
"usage: fasta_uncollapser [-c N] [-h] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n" \
"   [-c N]       = Assume input is a tabular file (not FASTA file),\n" \
"                  And the collapsed identifier (e.g. '1-1000') is on column N.\n" \
"   [-i INFILE]  = FASTA/Tabular input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Tabular output file. default is STDOUT.\n" \
"\n";

size_t collapsed_identifier_column = 0;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
		case 'c':
			if (optarg==NULL)
				errx(1, "[-c] parameter requires an argument value");
			collapsed_identifier_column = strtoul(optarg,NULL,10);
			if (collapsed_identifier_column<=0)
				errx(1,"Invalid column number (-c %s)", optarg);
			break;
		default:
			errx(1,"Internal error: unknown option '%c'",optc);
	}
	return 1;
}

void uncollapse_fasta_file()
{
	FASTX fastx;
	fastx_init_reader(&fastx, get_input_filename(),
		FASTA_ONLY, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	size_t seqid=1;
	while ( fastx_read_next_record(&fastx) ) {
		int count = get_reads_count(&fastx);
		for (int i=0;i<count;++i) {
			snprintf(fastx.name,sizeof(fastx.name),"%zu",seqid);
			++seqid;
			fastx_write_record(&fastx);
		}
	}

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu sequences (representing %zu reads)\n",
				num_input_sequences(&fastx), num_input_reads(&fastx));
		fprintf(get_report_file(), "Output: %zu sequences (representing %zu reads)\n",
				num_output_sequences(&fastx), num_output_reads(&fastx));
	}
}

//extract collapsed value, if any
size_t extract_collapsed_read_count(const std::string& text)
{
	char *endptr;

	std::string value = text;

	//Try 'N-NNNN' (as produced by fastx_collapser
	std::string::size_type dash_idx = value.find('-');

	// minus character found, and NOT the last character
	if (dash_idx != std::string::npos) {
		if ( (dash_idx+1)<text.length() )
			value = text.substr(dash_idx+1);
		else
			return 1; //last character is a minus: not a recognizable collapsed value, just return 1
	}

	size_t count = strtoul(value.c_str(), &endptr, 10);
	if ( count>0 && *endptr==0 ) //value converted successfuly, without surplus character?
		return count;

	return 1;
}

void uncollapse_tabular_file()
{
	ios::sync_with_stdio(false);
	size_t input_count=0;
	size_t output_count=0;

	string input_file(get_input_filename());
	if (input_file=="-")
		input_file="";
	string output_file(get_output_filename());
	if (output_file=="-")
		output_file="";

	InputStreamWrapper input(input_file);
	OutputStreamWrapper output(output_file);

	TextLineReader reader(input.stream());

	ostream &os = output.stream();
//	size_t seqid=1;
	while (reader.next_line()) {
		++input_count;
		vector<string> tokens;

		//Split input line into fields
		String_Tokenize(reader.line_string(), back_inserter(tokens), "\t");

		if (tokens.size()<collapsed_identifier_column) {
			cerr << "Input error in file '" << get_input_filename()
				<< "' line " << reader.line_number()
				<< ": got only " << tokens.size()
				<< " columns, but collapsed identifier column (-c) is "
				<< collapsed_identifier_column << endl;
			exit(1);
		}

		size_t count = extract_collapsed_read_count(tokens[collapsed_identifier_column-1]);
		output_count += count;

#if 0
		/* replace the collapse-identifier column with a sequencial counter */
		std::string prefix_fields;
		if (collapsed_identifier_column>1) {
			prefix_fields = join(tokens.begin(),tokens.begin()+(collapsed_identifier_column-1),"\t");
			prefix_fields += "\t";
		}

		std::string suffix_fields;
		if (collapsed_identifier_column<tokens.size()) {
			suffix_fields = "\t";
			suffix_fields += join(tokens.begin()+collapsed_identifier_column, tokens.end(),"\t");
		}

		for (size_t i=0;i<count;++i) {
			os << prefix_fields << seqid << suffix_fields << endl;
			++seqid;
		}
#else
		for (size_t i=0;i<count;++i)
			os << reader.line_string() << endl;
#endif
	}
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu lines (with collapsed sequence identifiers)\n", input_count);
		fprintf(get_report_file(), "Output: %zu lines\n", output_count);
	}
}

int main(int argc, char* argv[])
{
	ofstream output_file ;

	fastx_parse_cmdline(argc, argv, "c:", parse_program_args );

	if (collapsed_identifier_column==0)
		uncollapse_fasta_file();
	else
		uncollapse_tabular_file();

	return 0;
}
