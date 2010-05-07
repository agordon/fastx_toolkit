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
#include <assert.h>
#include <cstdlib>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <memory>
#include <limits.h>
#include <tr1/unordered_map>
#include <gtextutils/container_join.h>
#include <gtextutils/stream_wrapper.h>

#include "sequence_alignment.h"

#include <errno.h>
#include <err.h>

#include <config.h>

//#include "kmer_coding.h"
#include "adapter_hash.h"
#include "fastx_writer.h"
#include "sequence_writer.h"
#include "tile_detection.h"

#include "fastx.h"
#include "fastx_args.h"


#define MAX_ADAPTER_LEN 100

const char* usage=
"usage: fastx_clipper_v3 [-h] [-a ADAPTER] [-D] [-l N] [-n] [-d N] [-c] [-C] [-o] [-v] [-z] [-i INFILE] [-o OUTFILE]\n" \
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

using namespace std;
using namespace std::tr1;

enum {
	OPT_CLIPPED_SEQUENCES_FILE = CHAR_MAX+1,
	OPT_UNCLIPPED_SEQUENCES_FILE,
	OPT_ADAPTERS_FILE,
	OPT_TOO_SHORT_FILE,
	OPT_MISSING_START_THRESHOLD,
	OPT_MISSING_START_ACTION
};

struct option clipper_argv_options[] = {
	{"clipped",	1,	NULL,	OPT_CLIPPED_SEQUENCES_FILE},
	{"unclipped",	1,	NULL,	OPT_UNCLIPPED_SEQUENCES_FILE},
	{"adapters",	1,	NULL,	OPT_ADAPTERS_FILE},
	{"tooshort",    1,	NULL,	OPT_TOO_SHORT_FILE},
	{"missing-start-threshold", 1,	NULL,	OPT_MISSING_START_THRESHOLD},
	{"missing-start-action",    1,	NULL,	OPT_MISSING_START_ACTION}
};

//Parameters set by command line
string adapter;
static int verbose=0;
int debug_print_adapter_hash=0;
int debug_print_query_hash=0;
int debug_print_tile_detection=0;
int debug_print_region_detection=0;
int debug_adapter_verification=0;
int debug_show_clipping_result=0;
int debug_show_local_alignment_matrix=0;
int FASTQ_quality_offset=64;
static string input_filename="-";
static string output_filename="-";
int output_unclipped_sequences=1;
int output_clipped_sequences=1;

int max_tile_distance_difference = 1;
unsigned int min_length_after_clipping=15;
//unsigned int min_adapter_length_tile_match=15;
unsigned int min_adapter_covered_bases=10;
int min_local_alignment_score=7 ;


string clipped_sequences_filename;
string unclipped_sequences_filename;
string adapters_sequences_filename;
string tooshort_sequences_filename;
int k=5;

//How to handle this case?
//The end is matching to the adapter, the beggning is not
//query: TTTTG GGTCGNGAGAGCGGTTCAGCAGGAATGCCGA GACCGATCTCGTATGCCGTCATCTGCTTGGAAGAGAGGA
//adapt: AAAAA GATCGGAAGAGCGGTTCAGCAGGAATGCCGA TCT
//Three options:
//   Don't clip at all.
//   Clip the entire thing (
//   Clip from the detected/matched fragment only,
// skip_adapter_start is "5" in the above case.
unsigned int missing_start_threshold = 5 ;
enum MISSING_START_ACTION {
	MISSING_START_DONT_CLIP = 1,
	MISSING_START_CLIP_ALL,
	MISSING_START_CLIP_DETECTED
} ;
//MISSING_START_ACTION missing_start_action = MISSING_START_CLIP_DETECTED;
MISSING_START_ACTION missing_start_action = MISSING_START_CLIP_ALL;



//Global Variables

FASTX fastx;
SequenceHash AdapterHash;
AutoSequenceWriter clipped_sequences_writer;
AutoSequenceWriter unclipped_sequences_writer;
AutoSequenceWriter adapters_sequences_writer;
AutoSequenceWriter tooshort_sequences_writer;


size_t hamming_distance(const std::string &a, const std::string &b)
{
	size_t distance=0;
	size_t len = (a.length()>b.length()) ? b.length() : a.length();
	for (size_t i = 0; i< len; ++i){
		if ( a[i] != b[i] )
			distance++;
	}
	return distance;
}

//HalfLocalSequenceAlignment align;

void debug_print_query_kmers ( const string & query )
{
	int k_step = 1;
	for (unsigned int i=0;i<query.length()-k;i+=k_step) {
		const string query_kmer = query.substr(i,k);

		if (AdapterHash.kmer_exists(query_kmer)) {
			cerr << i << "\t" << query_kmer << "\t";
			cerr	<< "offset:\t" << AdapterHash.get_kmer_offset(query_kmer);
			cerr << endl;
		}
	}
}

/*
class QueryKmers {
public:
	unsigned int	query_offset;
	std::string	kmer;

	QueryKmers ( unsigned int _query_offset, const std::string& _kmer ) :
		query_offset(_query_offset), kmer(_kmer)
	{
	}
} ;
typedef std::vector<QueryKmers> query_kmers_vector;
*/

#if 0
struct tile_detection_results
{
	unsigned int query_start;
	unsigned int query_end;

	unsigned int adapter_start;
	unsigned int adapter_end;

	unsigned int num_consecutive_tiles;
	unsigned int num_bases_covered;

/*	unsigned int query_offset;
	unsigned int adapter_offset;
	unsigned int adapter_match_length;
	vector<bool> adapter_coverage;

	*/
};
#endif 
struct verification_results {
	unsigned int query_clip_offset ;
	unsigned int adapter_trim_offset ;
};

bool get_adapter_match_region ( const string & query, DetectedTileRegion /*output*/ &region )
{
	detect_tiles_vector_type detected_tiles;

	bool b = detect_tiles ( query, AdapterHash, detected_tiles );
	if (debug_print_tile_detection)
		cerr << "Tile-Detection: " << endl
			<< join(detected_tiles, "");
	if (!b) {
		if (debug_print_tile_detection)
			cerr << "Tile-Detection: not enough tiles detected. not clipping." << endl;
		return false;
	}

	detected_tile_region_vector_type detected_regions;

	b = detect_tile_regions ( detected_tiles, detected_regions ) ;
	if (debug_print_region_detection)
		cerr << "Tile-Region-Detection: " << endl
			<< join(detected_regions, "");
	if (!b) {
		if (debug_print_region_detection)
			cerr << "Tile-Region-Detection: not enough regions detected. not clipping." << endl;
		return false;
	}

	join_close_regions ( detected_regions ) ;
	if (debug_print_region_detection)
		cerr << "Tile-Region-Detection (after join): " << endl
			<< join(detected_regions, "");

	region = get_longest_region(detected_regions);
	if (debug_print_region_detection)
		cerr << "Tile-Region-Detection (longest region): " << endl
			<< region << endl;

	return true;
}

bool verify_matched_region ( const std::string & query,
				DetectedTileRegion /* in/out */ &  match_region )
{
	//Check the first tile of the aligned query and adapter.
	//If they don't match, force a local-alignment on the estimated alignment region.
	//Verify the estimated alignment region
	DetectedTileRegion aligned_region(match_region);
	unsigned int delta = min(aligned_region.adapter_start, aligned_region.query_start);
	aligned_region.query_start -= delta;
	aligned_region.adapter_start -= delta;
	const string adapter_fragment( adapter.substr(aligned_region.adapter_start, 5) );
	const string query_fragment ( query.substr(aligned_region.query_start,5));
//	int distance = hamming_distance ( query_fragment, adapter_fragment ) ;
//	if (distance>1) {
		if (debug_adapter_verification)
			cerr << "Adapter-Verification:" << endl
				<< "  verifing with local-alignment" << endl;

		if (aligned_region.query_start>0) aligned_region.query_start--;
		if (aligned_region.query_start>0) aligned_region.query_start--;
		if (aligned_region.adapter_start>0) aligned_region.adapter_start--;
		if (aligned_region.adapter_start>0) aligned_region.adapter_start--;

		const string adapter_estimated_region (
				adapter.substr(aligned_region.adapter_start,
					aligned_region.adapter_end - aligned_region.adapter_start) ) ;

		const string query_estimated_region (
				query.substr(aligned_region.query_start,
					aligned_region.query_end - aligned_region.query_start )) ;

		SequenceAlignment sa;
		sa.align(query_estimated_region, adapter_estimated_region);

		const SequenceAlignmentResults &r = sa.alignment_results();
		if (debug_adapter_verification) {
			r.print_score(cerr);
			r.print_aligned_coordinates(cerr);
			r.print_alignment(cerr);
		}

		if ( (size_t)(r.score*3/2) > r.match_count ) {
			//Use the local alignment match resutls
			match_region.query_start = aligned_region.query_start + r.query_start;
			match_region.query_end = aligned_region.query_start + r.query_end;

			match_region.adapter_start = aligned_region.adapter_start + r.target_start;
			match_region.adapter_end = aligned_region.adapter_start + r.target_end;

			if (debug_adapter_verification)
				cerr << " using local-alignment results. " << endl;
		} else {
			if (debug_adapter_verification)
				cerr << " NOT using local-alignment results. " << endl;
		}
//	}

	if ( (match_region.adapter_end - match_region.adapter_start) < min_adapter_covered_bases ) {
		// Only a small sub-sequence matches ?

		if ( (query.length() - match_region.query_end) < 2 ) {
			//At the very end of the query - that's OK
			if (debug_adapter_verification)
				cerr << "Adapter-Verification:" << endl
					<< " Allowing short match at the end of the query" << endl;
			return true;
		} else {
			if (debug_adapter_verification)
				cerr << "Adapter-Verification:" << endl
				       << "  adapter-covered-bases ("
					<<(match_region.adapter_end - match_region.adapter_start)
					<< ") smaller than min_adapter_covered_bases ("
					<< min_adapter_covered_bases
					<< "). not clipping."
					<< endl;
			return false;
		}
	}
	//The start of the adapter sequence is missing, but the middle is detected.
	//However, it is detected in the middle of the query-sequence.
	//What to do ?
	if (match_region.adapter_start > missing_start_threshold
			&&
		match_region.query_start > missing_start_threshold) {
		if (debug_adapter_verification)
			cerr << "Adapter-Verification:" << endl
				<< " found MISSING_START: ";

		switch (missing_start_action)
		{
		case MISSING_START_DONT_CLIP:
			if (debug_adapter_verification)
				cerr << " Action = Don't-Clip" << endl;
			return false;

		case MISSING_START_CLIP_ALL:
			{
			unsigned int delta = min(match_region.adapter_start, match_region.query_start);
			match_region.query_start -= delta;
			match_region.adapter_start -= delta;
			if (debug_adapter_verification)
				cerr << " Action = Clip-all. Adjusted match-region:" << endl
					<< " " << match_region << endl;
			return true;
			}

		case MISSING_START_CLIP_DETECTED:
			if (debug_adapter_verification)
				cerr << " Action = Clip-Only-Detected-Region" << endl;
			return true;

		default:
			errx(1,"Internal error: skip_start_action has invalid value (%d)", (int)missing_start_action);
		}
	}
	return true;
}

#if 0
bool verify_short_adapter_in_query(const std::string& query, tile_detection_results & tile_detection, verification_results & /*verification*/ )
{
	int expected_start_offset =
		(int)tile_detection.query_start - (int)tile_detection.adapter_start ;
	if (expected_start_offset<0)
		expected_start_offset=0;

	if (debug_adapter_verification)
		cerr << "verify short-check: trying local-alignment with last "
			<<  (query.length() - expected_start_offset)
			<< " nucleotides in the query string"
			<< endl;

	string query_fragment = query.substr(expected_start_offset);
	string adapter_fragment = adapter.substr(0,query_fragment.length()); //+2 = allow two deletions in the query, so add two nucleotides in the adapter

	/*align.align(query_fragment, adapter_fragment);

	if (debug_show_local_alignment_matrix)
		align.print_matrix();

	if (debug_adapter_verification)
		align.results().print();

	if (align.results().score >= min_local_alignment_score) {
		if (debug_adapter_verification)
			cerr << "verify short-check: local alignment returned valid score, adapter at offset "
				<< expected_start_offset + align.results().query_start
				<< endl;
		verification.query_clip_offset = expected_start_offset + align.results().query_start ;
		verification.adapter_trim_offset = align.results().target_start ;
		return true;
	}
*/
	return false;
}
#endif 

#if 0
bool verify_adapter_in_query(const std::string& query, tile_detection_results & tile_detection, verification_results &verification )
{
	int expected_start_offset =
		(int)tile_detection.query_start - (int)tile_detection.adapter_start ;
	if (expected_start_offset<0)
		expected_start_offset=0;

	verification.adapter_trim_offset = tile_detection.adapter_start ;
	verification.query_clip_offset = expected_start_offset;

	if (debug_adapter_verification)
		cerr << "Adapter-Verification: " << endl
			<< " query_start = " << tile_detection.query_start
			<< " query_end = " << tile_detection.query_end
			<< " adapter_start = " << tile_detection.adapter_start
			<< " adapter_end = " << tile_detection.adapter_end
			<< " Num_consecutive_tiles = " << tile_detection.num_consecutive_tiles
			<< " Expectd_start_offset = " << expected_start_offset
			<< endl;

	//Short adapter ?
	if ( (tile_detection.adapter_end - tile_detection.adapter_start) < min_adapter_length_tile_match) {

		size_t extend_adapter_start = tile_detection.adapter_start ;
		size_t extend_adapter_end =   tile_detection.adapter_end ;
		size_t extend_query_start =   tile_detection.query_start ;
		size_t extend_query_end   =   tile_detection.query_end ;

		//Try hard to find valid alignment -
		// extend the detected overlap region, and test it
		if ( extend_adapter_start > 0 ) {
			//extend it backwards
			size_t delta = (extend_adapter_start>extend_query_start)?extend_query_start:extend_adapter_start;
			extend_query_start -= delta;
			extend_adapter_start -= delta ;
		}

		if (extend_query_end < query.length()) {
			//extend it forwards, as much as possible (until the end of the query)
			size_t query_delta = query.length() - extend_query_end;
			size_t adapter_delta = adapter.length() - extend_adapter_end;
			size_t delta = (adapter_delta>query_delta)?query_delta:adapter_delta;
			extend_query_end += delta;
			extend_adapter_end += delta;
		}

		string extend_adapter_fragment = adapter.substr(extend_adapter_start, extend_adapter_end-extend_adapter_start);
		string extend_query_fragment = query.substr(extend_query_start, extend_query_end-extend_query_start);


		if (debug_adapter_verification)
			cerr << "Adapter-Verification: short-adater, after extension:" << endl
				<< "  extend_query_start = " << extend_query_start
				<< " extend_query_end = " << extend_query_end
				<< " extend_adapter_start = " << extend_adapter_start
				<< " extend_adapter_end = " << extend_adapter_end
				<< endl
				<< " extended_query = " << extend_query_fragment << endl
				<< " extended_adptr = " << extend_adapter_fragment << endl
				;

		size_t distance = hamming_distance ( extend_adapter_fragment, extend_query_fragment ) ;
		double score = (extend_adapter_fragment.length() - distance) / (double)extend_adapter_fragment.length();

		if (debug_adapter_verification)
			cerr << "Adapter-Verification: hamming-distance = " << distance
				<< " score = " << score
				<< endl;

/*		align.align(extend_query_fragment, extend_adapter_fragment);
		if (debug_show_local_alignment_matrix)
			align.print_matrix();
		if (debug_adapter_verification)
			align.results().print();*/

		if ( score > (3/4)) {
			verification.query_clip_offset = extend_query_start ;
			verification.adapter_trim_offset = extend_adapter_start ;
			return true;
		}

		//If it happens at the beginning of the adapter, and the end of the query sequence,
		if ( tile_detection.adapter_start < (size_t)k
				&&
			tile_detection.query_start >= (query.length() - min_adapter_length_tile_match - k )) {

			if (debug_adapter_verification)
				cerr << "Adapter-Verification: found SHORT-ADAPTER" << endl;

			//return verify_short_adapter_in_query ( query, tile_detection, verification ) ;
			return true;
		}

		if (debug_adapter_verification)
			cerr << "Adapter-Verification: detected adapter is too short (" <<(tile_detection.adapter_end - tile_detection.adapter_start)
				<< ") and does not match the end of the query. not clipping."
				<< endl;
		return false;
	}

	//The start of the adapter sequence is missing, but the middle is detected.
	//However, it is detected in the middle of the query-sequence.
	//What to do ?
	if (tile_detection.adapter_start > missing_start_threshold
			&&
		tile_detection.query_start > missing_start_threshold) {
		if (debug_adapter_verification)
			cerr << "Adapter-Verification: found MISSING_START" << endl;

		switch (missing_start_action)
		{
		case MISSING_START_DONT_CLIP:
			return false;

		case MISSING_START_CLIP_ALL:
			return true;

		case MISSING_START_CLIP_DETECTED:
			verification.query_clip_offset = tile_detection.query_start;
			verification.adapter_trim_offset = tile_detection.adapter_start;
			return true;

		default:
			errx(1,"Internal error: skip_start_action has invalid value (%d)", (int)missing_start_action);
		}
	}

	return true;
}
#endif

int parse_commandline(int argc, char* argv[])
{
	int c;
	int option_index;
	while ( (c=getopt_long(argc,argv,"vi:o:Q:D:a:k:M:", clipper_argv_options, &option_index)) != -1 ) {
		switch(c)
		{
			/* Standard Options */
		case 'i':
			input_filename = optarg;
			break;
		case 'o':
			output_filename = optarg;
			break;
		case 'Q':
			FASTQ_quality_offset = atoi(optarg);
			if (FASTQ_quality_offset<=0)
				errx(1,"Invalid FASTQ quality offset (%s). Value must be bigger than zero.", optarg);
			break;
		case 'v':
			verbose++;
			break;

			/* Clipper options */
		case 'a':
			adapter = optarg;
			break;
		case 'D':
			if (strcmp(optarg,"adapterhash")==0)
				debug_print_adapter_hash = 1;
			else
			if (strcmp(optarg,"queryhash")==0)
				debug_print_query_hash = 1 ;
			else
			if (strcmp(optarg,"tiledetection")==0)
				debug_print_tile_detection = 1;
			else
			if (strcmp(optarg,"regiondetection")==0)
				debug_print_region_detection = 1;
			else
			if (strcmp(optarg,"adapterverification")==0)
				debug_adapter_verification = 1;
			else
			if (strcmp(optarg,"showclip")==0)
				debug_show_clipping_result = 1 ;
			else
			if (strcmp(optarg,"localalignment")==0)
				debug_show_local_alignment_matrix = 1 ;
			else
				errx(1,"Unknown debug option '%s'", optarg);
			break;
		case 'k':
			k = atoi(optarg);
			if ( k<4 )
				errx(1,"Invalid k-mer length (-k %s). Must be a number larger than 3.", optarg);
			break;
		case 'M':
			min_adapter_covered_bases = atoi(optarg);
			break;

			/* Advanced Clipper Options */
		case OPT_CLIPPED_SEQUENCES_FILE:
			//clipped_sequences_writer = AutoSequenceWriter ( new SequenceWriter ( optarg ) ) ;
			clipped_sequences_filename = optarg;
			break;
		case OPT_UNCLIPPED_SEQUENCES_FILE:
			unclipped_sequences_filename = optarg;
			break;
		case OPT_ADAPTERS_FILE:
			adapters_sequences_filename = optarg;
			break;
		case OPT_TOO_SHORT_FILE:
			tooshort_sequences_filename = optarg;
			break;
		case OPT_MISSING_START_THRESHOLD:
			missing_start_threshold = atoi(optarg);
			break;
		case OPT_MISSING_START_ACTION:
			if (strcmp(optarg,"dontclip")==0)
				missing_start_action = MISSING_START_DONT_CLIP;
			else
			if (strcmp(optarg,"clipall")==0)
				missing_start_action = MISSING_START_CLIP_ALL;
			else
			if (strcmp(optarg,"clipdetected")==0)
				missing_start_action = MISSING_START_CLIP_DETECTED;
			else
				errx(1,"Invalid --missing-start-action argument '%s'", optarg);
			break;

		default:
			errx(1,"Unknown command line option -%c", optopt);
			break;
		}
	}
	return 1;
}

void verify_command_line()
{
	if (adapter.empty())
		errx(1,"missing adapter sequence (-a XXXXXXXXX)" );
}

void create_helper_files()
{
	SequenceWriter::TYPE t =
		fastx.write_fastq ? SequenceWriter::FASTQ : SequenceWriter::FASTA;

	if (!clipped_sequences_filename.empty())
		clipped_sequences_writer =
			AutoSequenceWriter (
				new SequenceWriter ( clipped_sequences_filename, t, "::CLIPPED" ) ) ;

	if (!unclipped_sequences_filename.empty())
		unclipped_sequences_writer =
			AutoSequenceWriter (
				new SequenceWriter ( unclipped_sequences_filename, t, "::UNCLIPPED" ) ) ;

	if (!adapters_sequences_filename.empty())
		adapters_sequences_writer =
			AutoSequenceWriter (
				new SequenceWriter ( adapters_sequences_filename, t, "::ADAPTERS" ) ) ;

	if (!tooshort_sequences_filename.empty())
		tooshort_sequences_writer =
			AutoSequenceWriter (
				new SequenceWriter ( tooshort_sequences_filename, t, "::TOOSHORT" ) ) ;
}


int main(int argc, char* argv[])
{
	int reads_count;
	bool clipped_normal ;
	bool clipped_too_short ;
	bool clipped_adapter_only;


	parse_commandline(argc, argv);
	verify_command_line();
	//assert_valid_kmer_size(k);
	//assert_kmer_coding_decoding();

	AdapterHash.set_sequence(adapter,k);
	if (debug_print_adapter_hash)
		AdapterHash.debug_print_hash(cerr);

	fastx_init_reader(&fastx, input_filename.c_str(),
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		FASTQ_quality_offset );

	fastx_init_writer(&fastx, output_filename.c_str(), OUTPUT_SAME_AS_INPUT, compress_output_flag());

	create_helper_files();

	while ( fastx_read_next_record(&fastx) ) {

		reads_count = get_reads_count(&fastx);

		string query ( fastx.nucleotides ) ;

		if (debug_print_query_hash)
			debug_print_query_kmers ( query ) ;

		clipped_normal = false;
		clipped_too_short = false;
		clipped_adapter_only = false;

		DetectedTileRegion match_region;

		string query_before_clipping ;
		string query_after_clipping;
		string clipped_adapter_sequence ;

		if (get_adapter_match_region(query, match_region)) {
			if (verify_matched_region(query, match_region)) {
				//Adapter found, need to clip,
				//but: where does the adapter fall?
				if (match_region.query_start==0)
					clipped_adapter_only = 1 ;
				else
				if (match_region.query_start<min_length_after_clipping)
					clipped_too_short = 1 ;
				else
					clipped_normal = 1;


				query_before_clipping = fastx.nucleotides;
				clipped_adapter_sequence = query_before_clipping.substr(match_region.query_start,match_region.query_end-match_region.query_start);
				query_after_clipping = query_before_clipping.substr(0,match_region.query_start);

				if (debug_show_clipping_result) {
					string query_pre_padding ;
					string adapter_pre_padding ;
					string detected_adapter_padding ( match_region.query_start, ' ');
					if ( match_region.adapter_start > match_region.query_start ) {
						query_pre_padding = string( (match_region.adapter_start - match_region.query_start), ' ');
					} else {
						adapter_pre_padding = string( (match_region.query_start - match_region.adapter_start), ' ');
					}
					cerr << "Clipping-Results: fastx.name " << endl;
					cerr << "Orig. Seq:  " << query_pre_padding << query_before_clipping << endl;

					cerr << "Found Adpr: " << query_pre_padding << detected_adapter_padding << clipped_adapter_sequence << endl ;
					cerr << "Orig Adpr:  " << adapter_pre_padding << adapter << endl;
//						adapter.substr(match_region.adapter_start, match_region.adapter_end-match_region.adapter_start) << endl ;
					cerr << "After Clip: " << query_pre_padding << query_after_clipping << endl;
				}
			}
		}

		bool not_clipped = ! ( clipped_too_short || clipped_adapter_only || clipped_normal );

		//Write Helper Files - unclipped sequences
		if (not_clipped && unclipped_sequences_writer.get() != NULL)
			unclipped_sequences_writer->write(fastx);
		if ((!not_clipped) && adapters_sequences_writer.get() != NULL)
			adapters_sequences_writer->write(fastx,match_region.query_start);
		if (clipped_too_short && tooshort_sequences_writer.get() != NULL)
			tooshort_sequences_writer->write(fastx, 0, match_region.query_start);
		if (clipped_normal && clipped_sequences_writer.get() != NULL)
			clipped_sequences_writer->write(fastx, 0, match_region.query_start);

		//Old-style clipping, for libfastx.c
		if (clipped_normal)
			fastx.nucleotides[match_region.query_start]=0;

		//Write Output File
		if ( ( not_clipped && output_unclipped_sequences )
		     ||
		     ( clipped_normal && output_clipped_sequences ) )
			fastx_write_record(&fastx);
	}
	//
	//Print verbose report
	/*
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
	*/

	return 0;
}
