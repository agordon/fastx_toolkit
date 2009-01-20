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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <err.h>

#include <config.h>

#include "fastx.h"
#include "fastx_args.h"

#define MAX_ADAPTER_LEN 100

const char* usage=
"usage: fastx_clipper [-h] [-a ADAPTER] [-s N] [-l N] [-n] [-d N] [-c] [-C] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
"\n" \
"version " VERSION "\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-a ADAPTER] = ADAPTER string. default is CTGTAGGCACCATCAATC.\n" \
"   [-s N]       = Max. number of mismatches allowed. default is 2.\n" \
"   [-l N]       = discard sequences shorter than N nucleotides. default is 5.\n" \
"   [-d N]       = Keep the adapter and N bases after it.\n" \
"                  (using '-d 0' is the same as not using '-d' at all. which is the default).\n" \
"   [-c]         = Discard non-clipped sequences (i.e. - keep only sequences which contained the adapter).\n" \
"   [-C]         = Discard clipped sequences (i.e. - keep only sequences which did not contained the adapter).\n" \
"   [-n]         = keep sequences with unknown (N) nucleotides. default is to discard such sequences.\n" \
"   [-v]         = Verbose - report number of sequences.\n" \
"                  If [-o] is specified,  report will be printed to STDOUT.\n" \
"                  If [-o] is not specified (and output goes to STDOUT),\n" \
"                  report will be printed to STDERR.\n" \
"   [-z]         = Compress output with GZIP.\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";

//Default adapter (Solexa's MODBAN) used at the Hannon Lab
char adapter[MAX_ADAPTER_LEN]="CTGTAGGCACCATCAATC";
int max_mismatches=2;
int min_length=5;
int discard_unknown_bases=1;
int keep_delta=0;
int discard_non_clipped=0;
int discard_clipped=0;


//Statistics for verbose report
unsigned int count_input=0 ;
unsigned int count_discarded_too_short=0; // see [-l N] option
unsigned int count_discarded_adapter_at_index_zero=0;  //empty sequences (after clipping)
unsigned int count_discarded_no_adapter_found=0; // see [-c] option
unsigned int count_discarded_adapter_found=0; // see [-C] option
unsigned int count_discarded_N=0; // see [-n]

FASTX fastx;

/*
	match_score
		calculates the matching score between two given strings.
	input - 
		sequence - sequence NULL-terminated string to match.
		start_offset - first character to match in <sequence>
		modban   - second NULL-terminated string to match. 
	output - 
		match score (number of identical characters in both string).

	Example:
		
		match_score("ABCDE",0,"AACD") returns 3, becuase:
			ABCDE
			AACD
			----
			MxMM = 3 (M = match, x = no match)

		match_score("ABCDE",2,"AACD") returns 0, because:
			ABCDE
			  AACD
			------
			XXXXX  = 0 (X = no match)
		(note that <start_offset>==2, so the first two characters in <sequence> are skipped).
		
		match_score("AATTGATCAGGACATAGGACAACTGTAGGCACCAT", 22, "CTGTAGGCACCATCAATC") returns 12, because:
			AATTGATCAGGACATAGGACAACTGTAGACACCAT
			                      CTGTAGGCACCATCAATC
			----------------------------------------
					      MMMMMMxMMMMMM      = 12
		(note that <start_offset>==26, so the first 26 characters in <sequence> are skipped).
*/
int match_score(const char *sequence, int start_offset, 
		const char* modban)
{
	int score = 0 ;
	
	while (start_offset-- && *sequence!=0)
		sequence++;
	
	while ( *sequence!=0 && *modban!=0 ) {
		score += (*sequence == *modban);
		sequence++;
		modban++;
	}
	
	return score;
}

/*
	best_match_offset - 
		given two strings, returns the best (if any)
		matching offset.

	input - 
		sequence, modban - two NULL terminated string to match.

	output - 
		-1  =  no match was found.
		>=0 = offset in <sequence> where <modban> matches with highest score.
*/
int best_match_offset(const char *sequence, const char* modban)
{

	int i, overlapping_characters, matching_characters, score=0;
	int best_index=-1, best_score=0;
	int sequence_length = strlen(sequence);
	int modban_length = strlen(modban);
		
	// printf("Seq=%s (%d chars)\n\n", sequence, sequence_length );
	//
	for (i=0;i<sequence_length;i++) {
		overlapping_characters = ((sequence_length - i)>modban_length) ? modban_length : (sequence_length- i);
		matching_characters = match_score(sequence, i, modban);

		//Too many mismatches?
		if ( overlapping_characters - matching_characters > max_mismatches) 
			continue;	


		// For tiny overlaps (when the number of overlapping characters
		// between the sequence and the modban is smaller than (max_mismatches*2)
		// We require a perfect match (= matching_characeters==overlapping_characters )
		if ( (overlapping_characters <= (max_mismatches*2) )
			&&
	             (matching_characters != overlapping_characters) )
		     	continue;

		score = matching_characters ;
		
		//DEBUG
		//printf("Ofs %2d, overlapping %2d, matching %2d, MisMatch %2d, score %2d\n",
		//	i, overlapping_characters, matching_characters, 
		//	overlapping_characters - matching_characters, score) ;
		
		if (score>best_score) {
			best_index = i;
			best_score = score;
		}
	}
	
	return best_index;
}

int parse_program_args(int optind, int optc, char* optarg)
{
	switch(optc) {
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
			//if (!valid_sequence_string(modban)) 
			//	errx(1,"Invalid adapter string (-a %s)", adapter);
			break ;
			
		case 's':
			if (optarg==NULL) 
				errx(1, "[-s] parameter requires an argument value");
			
			max_mismatches = strtoul(optarg, NULL, 10);
			if (max_mismatches<0 || max_mismatches>5) 
				errx(1,"Invalid number mismatches specified (-s %s). Allowed range: 0 to 5.", optarg);
			break;

		case 'l':
			if (optarg==NULL) 
				errx(1,"[-l] parameter requires an argument value");
			
			min_length = strtoul(optarg, NULL, 10);
			if (min_length<0) 
				errx(1,"Invalid minimum length specified (-l %s)", optarg);
			break;
			
		case 'n':
			discard_unknown_bases = 0 ;
			break;


	}
	return 1;
}

int parse_commandline(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "Ccd:a:s:l:n", parse_program_args);

	if (keep_delta>0) 
		keep_delta += strlen(adapter);
	return 1;
}


int main(int argc, char* argv[])
{
	int i;
	int reads_count;

	parse_commandline(argc, argv);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE);

	fastx_init_writer(&fastx, get_output_filename(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	while ( fastx_read_next_record(&fastx) ) {

		reads_count = get_reads_count(&fastx);

		count_input+= reads_count;

		//Find the best match with the adapter
		i = best_match_offset(fastx.nucleotides, adapter);
		if (i!=-1 && i>0) {
			i += keep_delta;
			//Just trim the string after this position
			fastx.nucleotides[i] = 0 ;
		}

		if (i==0) { // empty sequence ? (in which the adapter was found at index 0)
			count_discarded_adapter_at_index_zero += reads_count;
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
		
		//none of the above condition matched, so print this sequence.
		fastx_write_record(&fastx);
	}

	//
	//Print verbose report
	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Clipping Adapter: %s\n", adapter );
		fprintf(get_report_file(), "Max. mismatches: %d\n", max_mismatches ) ;
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
