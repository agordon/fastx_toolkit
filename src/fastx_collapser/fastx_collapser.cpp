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
#include <cassert>
#include <getopt.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <stdio.h>

#include <tr1/unordered_set>

#include "config.h"

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "libfastx/sequence.h"
#include "libfastx/fastx_file.h"
#include "libfastx/tab_file.h"
#include "libfastx/strings_buffer.h"
#include "libfastx/fastxse_commandline_parameters.h"
//#include <google/sparse_hash_map>

//#include "MurmurHash2_64.cpp"

using namespace std;

enum OPTIONS {
	OPT_SLEEP_WHEN_DONE = CHAR_MAX+1,
	OPT_OLD_METHOD
};

class FastxCollapserCommandLine : public FastxSE_commandline_parameters
{
public:
	size_t sleep;
	bool old_method;

	FastxCollapserCommandLine() :
		FastxSE_commandline_parameters(),
		sleep(0),
		old_method(false)
	{
		add_option_long("sleep", OPT_SLEEP_WHEN_DONE, required_argument);
		add_option_long("old-method", OPT_SLEEP_WHEN_DONE, no_argument);
	}

	void parameter_action(const int short_option_value, const std::string& optarg)
	{
		switch (short_option_value)
		{
		case OPT_SLEEP_WHEN_DONE:
			sleep = atoi(optarg.c_str());
			if (sleep==0)
				errx(1,"Error: invalid value for --sleep '%s'", optarg.c_str());
			break;

		case OPT_OLD_METHOD:
			old_method = true ;
			break;

		default:
			FastxSE_commandline_parameters::parameter_action(short_option_value,optarg);
			break;
		};
	}

	void print_help()
	{
	const char* usage=
"usage: fastx_collapser [-p] [-h] [-v] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h]         = This helpful help screen.\n" \
"   [-v]         = verbose: print short summary of input/output counts\n" \
"   [-p]         = Collapse sequences with identical prefix\n" \
"   [-i INFILE]  = FASTA/Q input file. default is STDIN.\n" \
"   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.\n" \
"\n";
		cout << usage ;
	}
};


#include <tr1/unordered_map>
std::tr1::unordered_map<string,size_t> collapsed_sequences;

typedef std::pair<string,size_t> sequence_count_pair;
typedef std::list< sequence_count_pair > sequence_count_list ;
sequence_count_list sorted_collapsed_sequences ;
bool flag_collapse_by_prefix=false;
bool flag_debug_collapse_by_prefix=false;

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

inline
bool sort_by_sequence ( const sequence_count_pair& sequence1, const sequence_count_pair& sequence2 )
{
	return sequence1.first < sequence2.first ;
}

inline
bool sort_by_abundance_count ( const sequence_count_pair& sequence1, const sequence_count_pair& sequence2 )
{
	return sequence1.second < sequence2.second ;
}


StringsBuffer buffer;

#define GCC_PACK __attribute__ ((packed))

class CollapsedSequenceData
{
private:
	unsigned int _multiplicity_count GCC_PACK;
	unsigned short _sequence_length GCC_PACK;
	char   _sequence[2] ;

	static size_t total_objects_allocated ;
	static size_t total_bytes_allocated ;

	CollapsedSequenceData();
	CollapsedSequenceData(const CollapsedSequenceData& other);

	CollapsedSequenceData(const Sequence& seq, bool add_quality=false )
	{
		_multiplicity_count = 1 ; //TODO - read real multiplicity count

		if (seq.nucleotides.length() > USHRT_MAX) {
			cerr << "Internal error: input sequence is longer than "
				<< USHRT_MAX << " nucleotides. Please recompile with larger variable type."
				<< endl;
			exit(1);
		}
		_sequence_length = seq.nucleotides.length();

		//Implicit assumption: enough memory was previously allocated for the entire
		//buffer!
		memcpy(_sequence, seq.nucleotides.data(), _sequence_length);
		_sequence[_sequence_length] = 0;

		if (add_quality) {
			assert(seq.quality_cached_line.length() == seq.nucleotides.length());

			char* quality_ptr = _sequence + _sequence_length +1 ;
			memcpy(quality_ptr, seq.quality_cached_line.data(), _sequence_length);
			quality_ptr[_sequence_length] = 0;
		}
	}

public:
	// The only way to allocate + create a CollapsedSequenceData
	static CollapsedSequenceData* create( StringsBuffer &buffer, const Sequence& seq, bool add_quality=false )
	{
		size_t needed_size = sizeof(CollapsedSequenceData)-1 +
				     seq.nucleotides.length() + 1;
		if (add_quality)
			needed_size += seq.quality_cached_line.length()+1;

		++total_objects_allocated;
		total_bytes_allocated += needed_size ;

		void *ptr = buffer.allocate_buffer0(needed_size);
		CollapsedSequenceData* csd = new (ptr) CollapsedSequenceData(seq, add_quality);
		return csd;
	}

	unsigned int multiplicity_count() const { return _multiplicity_count; }

	void add_multiplicity_count(const unsigned int count) { _multiplicity_count += count ; }

	size_t sequence_length() const { return _sequence_length ; }
	const char* sequence() const { return _sequence; }

	const char* quality() const { return (char*)(_sequence + _sequence_length + 1 ); }

	static size_t objects_count() { return total_objects_allocated; }
	static size_t bytes_count() { return total_bytes_allocated; }
};

struct CollapsedSequenceData_Hash
{
	size_t operator()(const CollapsedSequenceData* csd) const
	{
		return std::tr1::hash<const char*>()(csd->sequence());
	}
};

struct CollapsedSequenceData_EqualTo : public equal_to<CollapsedSequenceData*>
{
	bool operator()(const CollapsedSequenceData* a, const CollapsedSequenceData* b) const
	{
		if (a->sequence_length() != b->sequence_length())
			return false;
		return strncmp(a->sequence(), b->sequence(), a->sequence_length());
	}

	bool operator()(const CollapsedSequenceData* a, const std::string& b) const
	{
		if (a->sequence_length() != b.length())
			return false;
		return strncmp(a->sequence(), b.c_str(), a->sequence_length());
	}
};

struct eqstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	}
};

struct hash_char_ptr
{
	size_t operator()(const char* str) const
	{
		return std::tr1::hash<const char*>()(str);
	}
};


size_t CollapsedSequenceData::total_objects_allocated = 0 ;
size_t CollapsedSequenceData::total_bytes_allocated = 0;

typedef std::tr1::unordered_map<const char*,CollapsedSequenceData* ,hash_char_ptr, eqstr> CollapsedSequencesHash ;
//typedef google::sparse_hash_map<const char*,CollapsedSequenceData* ,hash_char_ptr, eqstr> CollapsedSequencesHash ;
typedef std::pair<const char*, CollapsedSequenceData*> CSD_PAIR;
CollapsedSequencesHash hash;

/*
struct equal_collapsed_sequence_data
{
	bool operator()(const CollapsedSequenceData* s1, const CollapsedSequenceData* s2) const
	{
		return (s1 == s2) ||
			(
			 (s1 && s2)
		          &&
			 (s1->seq_length == s2->seq_length)
			  &&
			 (strcmp(s1->sequence, s2->sequence) == 0)
			);
	}
};

struct murmur_hash_collapsed_sequence_data
{
	size_t operator()(const CollapsedSequenceData* s) const
	{
		return MurmurHash64A(s->sequence, s->seq_length, 0);
	}
};


typedef std::tr1::unordered_map<const char*, CollapsedSequenceData*> CollapsedSequencesMap;
CollapsedSequencesMap collapsed_sequences_map;

void add_collapsed_sequence(const Sequence &s)
{
	CollapsedSequencesMap::iterator it = collapsed_sequences_map.find(s.nucleotides.c_str());
	if (it == collapsed_sequences_map.end()) {
		//create new
		size_t size = sizeof(CollapsedSequenceData)-1 +
				s.nucleotides.length()+1 ;

		CollapsedSequenceData *pData = (CollapsedSequenceData*)buffer.allocate_buffer0(size);
		pData->multiplicity_count = 1 ;
		strcpy(pData->sequence, s.nucleotides.c_str());

		collapsed_sequences_map.insert(std::pair<const char*, CollapsedSequenceData*>(pData->sequence, pData));
	} else {
		//increment multiplicity-count
		it->second->multiplicity_count += 1 ;
	}
}
*/

#if 0
int main(int argc, char* argv[])
{
	ofstream output_file ;

	fastx_parse_cmdline(argc, argv, "p", parse_program_args );

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


	if (flag_collapse_by_prefix) {
		sorted_collapsed_sequences.sort ( sort_by_sequence ) ;

		sequence_count_list::iterator curr = sorted_collapsed_sequences.begin();
		while ( curr != sorted_collapsed_sequences.end() ) {
			const string& curr_seq = curr->first;

			sequence_count_list::iterator next = curr;
			next++;
			if (next==sorted_collapsed_sequences.end())
				break;

			const string& next_seq = next->first ;

			if (flag_debug_collapse_by_prefix)
			cerr << "--Checking:" << endl
				<< "curr: " << curr_seq << "\t" << curr->second << endl
				<< "next: " << next_seq << "\t" << next->second << endl;
			if ( (next_seq.length() >= curr_seq.length()
			      &&
			      next_seq.substr(0,curr_seq.length()) == curr_seq)
			      ||
			      (curr_seq.length() > next_seq.length()
			       &&
			       curr_seq.substr(0,next_seq.length()) == next_seq)
			      ) {
					if (flag_debug_collapse_by_prefix)
						cerr << "Collapsing by prefix, increasing next's count from " << next->second << " to " << (next->second + curr->second) << endl;

					next->second += curr->second;
					if (curr_seq.length() > next_seq.length())
						next->first = curr->first;
					curr = sorted_collapsed_sequences.erase(curr);
			}
			else
				curr = next ;
		}
	}

	sorted_collapsed_sequences.sort ( sort_by_abundance_count ) ;

	PrintCollapsedSequence stats =  for_each ( sorted_collapsed_sequences.rbegin(), 
			sorted_collapsed_sequences.rend(), PrintCollapsedSequence(real_output) ) ;

	/* This (in)sanity check prevents collapsing an already-collapsed FASTA file, so skip it for now */
	//if (stats.total_reads != num_input_reads(&fastx))
	//	errx(1,"Internal error: stats.total_reads (%zu) != num_input_reads(&fastx) (%zu).\n", 
	//		stats.total_reads, num_input_reads(&fastx) ); 
	//

	if ( verbose_flag() ) {
		fprintf(get_report_file(), "Input: %zu sequences (representing %zu reads)\n",
				num_input_sequences(&fastx), num_input_reads(&fastx));
		fprintf(get_report_file(), "Output: %zu sequences (representing %zu reads)\n",
				stats.counter, stats.total_reads);
	}
	return 0;
}
#endif

int main(int argc, char* argv[] )
{
//	cout << "sizeof(CollapsedSequenceData) = " << sizeof(CollapsedSequenceData) << endl;
//	exit(1);

	FastxCollapserCommandLine p;
	p.parse_command_line(argc,argv);

	Sequence seq;
	ISequenceReader *reader = p.reader().get();
	//ISequenceWriter *writer = p.writer().get();

	while (reader->read_next_sequence(seq)) {
		CollapsedSequencesHash::iterator it = hash.find(seq.nucleotides.c_str());
		if (it == hash.end()) {
			CollapsedSequenceData *ptr = CollapsedSequenceData::create(buffer, seq, false);
			CSD_PAIR p = CSD_PAIR( ptr->sequence(), ptr ) ;
			hash.insert(p);
		} else {
			CollapsedSequenceData* csd = it->second;
			csd->add_multiplicity_count(1);
		}
	}

	/*
	size_t i = 1 ;
	for ( sequence_count_list::reverse_iterator it  = sorted_collapsed_sequences.rbegin();
			it != sorted_collapsed_sequences.rend(); ++it ) {
		seq.nucleotides = it->first;
		stringstream ss;
		ss << i << "-" << it->second ;
		seq.id = ss.str();
		++i;

		writer->write_sequence(seq);
	}
	*/
	if (p.sleep) {
		cerr << CollapsedSequenceData::objects_count() << " objects, "
			<< CollapsedSequenceData::bytes_count() << " bytes." << endl;
		cerr << "Done! going to sleep for " << p.sleep << " seconds" << endl;
		sleep(p.sleep);
	}
	return 1;
}
