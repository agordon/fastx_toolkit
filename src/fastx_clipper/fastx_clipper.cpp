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
#include <algorithm>
#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>

#include "sequence_alignment.h"

#include <errno.h>
#include <err.h>

#include <config.h>

#include "fastx.h"
#include "fastx_args.h"


#define MAX_ADAPTER_LEN 100

const char* usage=
"usage: fastx_clipper [-h] [-a ADAPTER] [-D] [-l N] [-n] [-d N] [-c] [-C] [-o] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-a ADAPTER] = ADAPTER string. default is CCTTAAGG (dummy adapter).\n" \
"   [-l N]       = discard sequences shorter than N nucleotides. default is 5.\n" \
"   [-d N]       = Keep the adapter and N bases after it.\n" \
"                  (using '-d 0' is the same as not using '-d' at all. which is the default).\n" \
"   [-c]         = Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).\n" \
"   [-C]         = Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).\n" \
"   [-k]         = Report Adapter-Only sequences.\n" \
"   [-n]         = keep sequences with unknown (N) nucleotides. default is to discard such sequences.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-D]	 = DEBUG output.\n" \
"   [-M N]       = require minimum adapter alignment length of N.\n" \
"                  If less than N nucleotides aligned with the adapter - don't clip it." \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

//Default adapter - Dummy sequence
char adapter[MAX_ADAPTER_LEN]="CCTTAAGG";
unsigned int min_length=5;
int discard_unknown_bases=1;
int keep_delta=0;
int discard_non_clipped=0;
int discard_clipped=0;
int show_adapter_only=0;
int debug = 0 ;
int minimum_adapter_length = 0;


//Statistics for verbose report
unsigned int count_input=0 ;
unsigned int count_discarded_too_short=0; // see [-l N] option
unsigned int count_discarded_adapter_at_index_zero=0;  //empty sequences (after clipping)
unsigned int count_discarded_no_adapter_found=0; // see [-c] option
unsigned int count_discarded_adapter_found=0; // see [-C] option
unsigned int count_discarded_N=0; // see [-n]

FASTX fastx;
HalfLocalSequenceAlignment align;

int parse_program_args(int __attribute__((unused)) optind, int optc, char* optarg)
{
	switch(optc) {
		case 'M':
			if (optarg==NULL) 
				errx(1, "[-M] parameter requires an argument value");
			minimum_adapter_length = atoi(optarg);
			if (minimum_adapter_length<=0) 
				errx(1,"Invalid minimum adapter length (-M %s)", optarg);
			break;

		case 'k':
			show_adapter_only=1;
			break;

		case 'D':
			debug++;
			break ;

		case 'c':
			discard_non_clipped = 1;
			break;

		case 'C':
			discard_clipped = 1 ;
			break ;
		case 'd':
			if (optarg==NULL) 
				errx(1, "[-d] parameter requires an argument value");
			keep_delta = strtoul(optarg,NULL,10);
			if (keep_delta<0) 
				errx(1,"Invalid number bases to keep (-d %s)", optarg);
			break;
		case 'a':
			strncpy(adapter,optarg,sizeof(adapter)-1);
			//TODO:
			//if (!valid_sequence_string(adapter)) 
			//	errx(1,"Invalid adapter string (-a %s)", adapter);
			break ;
			
		case 'l':
			if (optarg==NULL) 
				errx(1,"[-l] parameter requires an argument value");
			
			min_length = strtoul(optarg, NULL, 10);
			break;
			
		case 'n':
			discard_unknown_bases = 0 ;
			break;

		default:
			errx(1,"Unknown argument (%c)", optc ) ;

	}
	return 1;
}

int parse_commandline(int argc, char* argv[])
{

	fastx_parse_cmdline(argc, argv, "M:kDCcd:a:s:l:n", parse_program_args);

	if (keep_delta>0) 
		keep_delta += strlen(adapter);
	return 1;
}

int adapter_cutoff_index ( const SequenceAlignmentResults& alignment_results ) __attribute__ ((const));
int adapter_cutoff_index ( const SequenceAlignmentResults& alignment_results )
{
	#if 0
	int mismatches = alignment_results.mismatches ;
	
	//The adapter(=target) is expected to align from the first base.
	//If the start is not zero (=not aligned from first base),
	//count each skipped base as a mismatch
	mismatches += alignment_results.target_start ;

	//The adapter is expected to align up to the end
	//of the adapter(=target), or the end of the query.
	//If it doesn't, count the un-aligned bases as mismatches
	int missing_from_query_end = (alignment_results.query_size - alignment_results.query_end-1);
	int missing_from_target_end = (alignment_results.target_size - alignment_results.target_end-1);

	int missing_from_end = std::min(missing_from_query_end, missing_from_target_end);
	
	mismatches += missing_from_end ;
	

	 
	std::cout << "Missing from start = " << alignment_results.target_start
		  << " Missing from end = " << missing_from_end
		  << " mismatches = " << mismatches 
		  << std::endl;
	
	if (mismatches > max_mismatches)
		return -1;

	return alignment_results.query_start;
	#endif

	int alignment_size = alignment_results.neutral_matches +
			     alignment_results.matches + 
			     alignment_results.mismatches +
			     alignment_results.gaps ;

	//No alignment at all?
	if (alignment_size==0)
		return -1;

	if (minimum_adapter_length>0 && alignment_size<minimum_adapter_length)
		return -1;

	//Any good alignment at the end of the query
	//(even only a single nucleotide)
	//Example:
	//  The adapter starts with CTGTAG, The Query ends with CT - it's a match.
	if ( alignment_results.query_end == alignment_results.query_size-1
	     &&
	     alignment_results.mismatches == 0 ) {
	     	//printf("--1\n");
		return alignment_results.query_start ;
	}

	if ( alignment_size > 5
	     &&
	     alignment_results.target_start == 0
	     &&
	     (alignment_results.matches * 100 / alignment_size ) >= 75 ) {
	     	//printf("--2\n");
		return alignment_results.query_start ;
	}

	if ( alignment_size > 11 
	     &&
	     (alignment_results.matches * 100 / alignment_size ) >= 80 ) {
	     	//printf("--2\n");
		return alignment_results.query_start ;
	}

	//
	//Be very lenient regarding alignments at the end of the query sequence
	if ( alignment_results.query_end >= alignment_results.query_size-2
	     &&
	     alignment_size <= 5 && alignment_results.matches >= 3) {
			//printf("--3\n");
			return alignment_results.query_start ;
		}

	return -1;
}


int main(int argc, char* argv[])
{
	int i;
	int reads_count;

	parse_commandline(argc, argv);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		reads_count = get_reads_count(&fastx);
		
		#if 0
		std::string query = std::string(fastx.nucleotides) + std::string( strlen(adapter), 'N' ); 
		std::string target= std::string( strlen(fastx.nucleotides), 'N' ) + std::string(adapter);
		#else
		std::string query = std::string(fastx.nucleotides) ;
		std::string target= std::string(adapter);
		#endif
		
		
		align.align( query, target ) ;

		if (debug>1) 
			align.print_matrix();
		if (debug>0)
			align.results().print();
		
		count_input+= reads_count;

		//Find the best match with the adapter
		i = adapter_cutoff_index ( align.results() ) ;
		
		if (i!=-1 && i>0) {
			i += keep_delta;
			//Just trim the string after this position
			fastx.nucleotides[i] = 0 ;
		}

		if (i==0) { // empty sequence ? (in which the adapter was found at index 0)
			count_discarded_adapter_at_index_zero += reads_count;
			
			if (show_adapter_only)
				fastx_write_record(&fastx);
			continue;
		}

		if (strlen(fastx.nucleotides) < min_length) { // too-short sequence ?
			count_discarded_too_short += reads_count;
			continue;
		}

		if ( (i==-1) && discard_non_clipped ) { // adapter not found (i.e. sequence was not clipped) ?
			count_discarded_no_adapter_found += reads_count;
			continue ;
		}

		if ( (i>0) && discard_clipped ) { // adapter found, and user requested to keep only non-clipped sequences 
			count_discarded_adapter_found += reads_count;
			continue;
		}

		if ( (discard_unknown_bases && strchr(fastx.nucleotides,'N')!=NULL ) ) { // contains unknown bases (after clipping) ?
			count_discarded_N += reads_count;
			continue;
		}
		
		if (!show_adapter_only)  {
			//none of the above condition matched, so print this sequence.
			fastx_write_record(&fastx);
		}
	}

	//
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Clipping Adapter: %s\n", adapter );
		fprintf(get_report_file(), "Min. Length: %d\n", min_length) ;

		if (discard_clipped)
			fprintf(get_report_file(), "Clipped reads - discarded.\n"  ) ;
		if (discard_non_clipped)
			fprintf(get_report_file(), "Non-Clipped reads - discarded.\n"  ) ;

		
		fprintf(get_report_file(), "Input: %u reads.\n", count_input ) ;
		fprintf(get_report_file(), "Output: %u reads.\n", 
			count_input - count_discarded_too_short - count_discarded_no_adapter_found - count_discarded_adapter_found -
			count_discarded_N - count_discarded_adapter_at_index_zero ) ;

		fprintf(get_report_file(), "discarded %u too-short reads.\n", count_discarded_too_short ) ;
		fprintf(get_report_file(), "discarded %u adapter-only reads.\n", count_discarded_adapter_at_index_zero );
		if (discard_non_clipped)
			fprintf(get_report_file(), "discarded %u non-clipped reads.\n", count_discarded_no_adapter_found );
		if (discard_clipped)
			fprintf(get_report_file(), "discarded %u clipped reads.\n", count_discarded_adapter_found );
		if (discard_unknown_bases)
			fprintf(get_report_file(), "discarded %u N reads.\n", count_discarded_N );
	}

	return 0;
}
