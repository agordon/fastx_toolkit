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
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <map>
#include <list>
#include <stdio.h>

#include "config.h"

#include "fastx.h"
#include "fastx_args.h"

using namespace std;

const char* usage=
"usage: fastx_collapser [-h] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

FASTX fastx;
#include <tr1/unordered_map>
std::tr1::unordered_map<string,size_t> collapsed_sequences;
std::list< pair<string,size_t> > sorted_collapsed_sequences ;

struct PrintCollapsedSequence
{
	size_t counter;
	size_t total_reads ;

	ostream &output ;
	PrintCollapsedSequence( ostream& _output ) : 
		counter(0), 
		total_reads(0),
		output(_output) {}

	void operator() ( const std::pair<string, int> & sequence )
	{
		counter++;
		total_reads += sequence.second ;
		output << ">" << counter << "-" << sequence.second << endl << sequence.first << endl ;
	}
};

bool sort_by_abundance_count ( const pair<string, size_t> & sequence1, const pair<string, size_t>& sequence2 )
{
	return sequence1.second < sequence2.second ;
}

int main(int argc, char* argv[])
{
	ofstream output_file ;

	fastx_parse_cmdline(argc, argv, "", NULL );

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	bool use_stdout = true;
	if ( strcmp(get_output_filename(), "-")!=0 ) {
		use_stdout = false;
		output_file.open(get_output_filename());
		if (!output_file) 
			errx(1,"Failed to create output file (%s)", get_output_filename() );
	}
	ostream& real_output = (use_stdout) ? cout : output_file ;

	while ( fastx_read_next_record(&fastx) ) {
		collapsed_sequences[string(fastx.nucleotides)]+= get_reads_count(&fastx);
	}
	
	copy ( collapsed_sequences.begin(), collapsed_sequences.end(), 
		back_inserter(sorted_collapsed_sequences) ) ;

	sorted_collapsed_sequences.sort ( sort_by_abundance_count ) ;

	PrintCollapsedSequence stats =  for_each ( sorted_collapsed_sequences.rbegin(), 
			sorted_collapsed_sequences.rend(), PrintCollapsedSequence(real_output) ) ;

	/* This (in)sanity check prevents collapsing an already-collapsed FASTA file, so skip it for now */
	/*
	if (stats.total_reads != num_input_reads(&fastx))
		errx(1,"Internal error: stats.total_reads (%zu) != num_input_reads(&fastx) (%zu).\n", 
			stats.total_reads, num_input_reads(&fastx) ); 
	*/

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu sequences (representing %zu reads)\n",
				num_input_sequences(&fastx), num_input_reads(&fastx));
		fprintf(get_report_file(), "Output: %zu sequences (representing %zu reads)\n",
				stats.counter, stats.total_reads);
	}
	return 0;
}
