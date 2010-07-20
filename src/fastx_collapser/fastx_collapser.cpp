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

#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include "config.h"

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "libfastx/sequence.h"
#include "libfastx/fastx_file.h"
#include "libfastx/tab_file.h"
#include "libfastx/strings_buffer.h"
#include "libfastx/fastxse_commandline_parameters.h"
#include "libfastx/MurmurHash2.h"

using namespace std;

enum {
	OPT_NO_SORT = CHAR_MAX+1,
	OPT_SORT_BY_SEQUENCE,
	OPT_SORT_BY_MULTIPLICITY
};

enum SORT_TYPE {
	SORT_BY_SEQUENCE,
	SORT_BY_MULTIPLICITY
} ;

class FastxCollapserCommandLine : public FastxSE_commandline_parameters
{
public:
	bool sort ;
	SORT_TYPE sort_type;

	FastxCollapserCommandLine() :
		FastxSE_commandline_parameters(),
		sort(true),
		sort_type(SORT_BY_MULTIPLICITY)
	{
		add_option_long("nosort", OPT_NO_SORT, no_argument);
		add_option_long("sortseq", OPT_SORT_BY_SEQUENCE, no_argument);
		add_option_long("sortmul", OPT_SORT_BY_MULTIPLICITY, no_argument);
	}

	void parameter_action(const int short_option_value, const std::string& optarg)
	{
		switch(short_option_value)
		{
		case OPT_NO_SORT:
			sort = false;
			break;

		case OPT_SORT_BY_SEQUENCE:
			sort = true;
			sort_type = SORT_BY_SEQUENCE;
			break;

		case OPT_SORT_BY_MULTIPLICITY:
			sort = true;
			sort_type = SORT_BY_MULTIPLICITY;
			break;

		default:
			FastxSE_commandline_parameters::parameter_action(short_option_value,optarg);
		}
	}

	void print_help()
	{
	const char* usage=
"usage: fastx_collapser [OPTIONS] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"options:\n" \
"   -i INFILE  = FASTA/Q input file. default is STDIN.\n" \
"   -o OUTFILE = FASTA/Q output file. default is STDOUT.\n" \
"   -h         = This helpful help screen.\n" \
"   -v         = verbose: print short summary of input/output counts.\n" \
"   --nosort   = Do no sort the output (saves some processing time).\n" \
"   --sortmul  = Sort the sequences by descending multiplicity counts (default).\n" \
"   --sortseq  = Sort the sequences by nucleotides string.\n" \
"   --pipein\n" \
"   --tabin    = Input is tabular file, not a FASTA/Q file.\n" \
"   --pipeout\n" \
"   --tabout   = Output tabular file, not a FASTA/Q file.\n" \
"\n";
		cout << usage ;
	}
};

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
		_multiplicity_count = seq.get_multiplicity_count();

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

	char* ptr_quality() { return (char*)(_sequence + _sequence_length + 1 ); }

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

	void merge_quality_scores ( const Sequence& s )
	{
		if (s.quality_cached_line.length() != sequence_length()) {
			cerr << "Internal error: trying to merge quality scores of different lengths (sequence()=" << sequence() << " sequence_length()=" << sequence_length() << " new_seq.nucleotides()=" << s.nucleotides << " new_seq.qualities.lengh()=" << s.quality_cached_line.length() << ")" << endl;
			exit(1);
		}

		char *qual_ptr = ptr_quality();
		for (size_t i=0;i<sequence_length();++i) {
			const char a = qual_ptr[i];
			const char b = s.quality_cached_line[i];
			if (b>a)
				qual_ptr[i]=b;
		}
	}

	size_t sequence_length() const { return _sequence_length ; }
	const char* sequence() const { return _sequence; }

	const char* quality() const { return (char*)(_sequence + _sequence_length + 1 ); }

	static size_t objects_count() { return total_objects_allocated; }
	static size_t bytes_count() { return total_bytes_allocated; }

};

struct CollapsedSequenceData_CompareSequence
{
	bool operator()(const CollapsedSequenceData* a, const CollapsedSequenceData* b) const
	{
		return strcmp(a->sequence(), b->sequence())<0;
	}
};

struct CollapsedSequenceData_CompareMultiplicity
{
	bool operator()(const CollapsedSequenceData* a, const CollapsedSequenceData* b) const
	{
		if (a->multiplicity_count()==b->multiplicity_count())
			return strcmp(a->sequence(), b->sequence())<0;

		return a->multiplicity_count() > b->multiplicity_count();
	}
};

struct eqstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	}
};

struct MurmurHash2_CSD
{
	size_t operator()(const CollapsedSequenceData* csd) const
	{
		return MurmurHash64A(csd->sequence(), csd->sequence_length(),0);
	}
};

struct EqualTo_CSD
{
	bool operator()(const CollapsedSequenceData* a, const CollapsedSequenceData* b) const
	{
		return ( (a==b)
			 ||
			 (
			  (a->sequence_length() == b->sequence_length())
			   &&
			  (strncmp(a->sequence(), b->sequence(), a->sequence_length())==0)
			 )
			);
	}
};

size_t CollapsedSequenceData::total_objects_allocated = 0 ;
size_t CollapsedSequenceData::total_bytes_allocated = 0;

#define USE_MAP

#ifdef USE_MAP
typedef std::tr1::unordered_map<const char*,CollapsedSequenceData* ,MurmurHash2_char_ptr, eqstr> CollapsedSequencesHash ;
#endif

#ifdef USE_SET
typedef std::tr1::unordered_set<CollapsedSequenceData*, MurmurHash2_CSD, EqualTo_CSD> CollapsedSequencesHash;
#endif

CollapsedSequencesHash hash;


int main(int argc, char* argv[] )
{
	FastxCollapserCommandLine p;
	p.parse_command_line(argc,argv);

	Sequence seq;
	ISequenceReader *reader = p.reader().get();
	ISequenceWriter *writer = p.writer().get();

	SequenceCounter input_counter;
	bool have_quality_scores=false;

	//Load data, collapse by sequence
#ifdef USE_MAP
	while (reader->read_next_sequence(seq)) {
		input_counter.add(seq);
		CollapsedSequencesHash::iterator it = hash.find(seq.nucleotides.c_str());
		if (it == hash.end()) {
			have_quality_scores = seq.have_quality_scores;
			CollapsedSequenceData *ptr = CollapsedSequenceData::create(buffer, seq, have_quality_scores);
			hash.insert( make_pair ( ptr->sequence(), ptr) );
		} else {
			CollapsedSequenceData* csd = it->second;
			csd->add_multiplicity_count(seq.get_multiplicity_count());
			if (have_quality_scores)
				csd->merge_quality_scores(seq);
		}
	}
#endif

#ifdef USE_SET
	while (reader->read_next_sequence(seq)) {
		input_counter.add(seq);
		have_quality_scores = seq.have_quality_scores;
		CollapsedSequenceData *ptr = CollapsedSequenceData::create(buffer, seq, have_quality_scores);
		CollapsedSequencesHash::iterator it = hash.find(ptr);
		if (it == hash.end()) {
			hash.insert( ptr );
		} else {
			buffer.rollback();
			CollapsedSequenceData* csd = *it;
			csd->add_multiplicity_count(seq.get_multiplicity_count());
			if (have_quality_scores)
				csd->merge_quality_scores(seq);
		}
	}
#endif

	//Sort by Multiplicity
	typedef std::vector<CollapsedSequenceData*> CollapsedSequenceData_Vector;
	CollapsedSequenceData_Vector v;
	v.resize(hash.size());
	CollapsedSequenceData_Vector::iterator it_v = v.begin();
	CollapsedSequencesHash::const_iterator it_h = hash.begin();
	for ( ; it_h != hash.end(); ++it_v, ++it_h)
#ifdef USE_SET
		*it_v = *it_h;
#endif
#ifdef USE_MAP
		*it_v = it_h->second;
#endif

	if (p.sort) {
		switch(p.sort_type)
		{
		case SORT_BY_SEQUENCE:
			sort(v.begin(), v.end(), CollapsedSequenceData_CompareSequence() );
			break;

		case SORT_BY_MULTIPLICITY:
			sort(v.begin(), v.end(), CollapsedSequenceData_CompareMultiplicity() ) ;
			break;
		default:
			cerr << "Internal Error: invalid sort-type (" << p.sort_type << ")" << endl;
			exit(1);
			break;
		}
	}

	//Write sequences
	SequenceCounter output_counter;
	size_t count = 1;
	for (it_v = v.begin(); it_v != v.end(); ++it_v) {
		const CollapsedSequenceData* ptr = *it_v;

		seq.clear();
		stringstream ss;
		ss << count << "-" << ptr->multiplicity_count();
		seq.id = ss.str();
		seq.nucleotides = ptr->sequence();

		if (have_quality_scores) {
			seq.quality_cached_line = ptr->quality();
			seq.id2 = seq.id;
			seq.have_quality_scores = true;
		}

		output_counter.add(seq);
		writer->write_sequence(seq);
		++count;
	}

	if (p.verbose()) {
		p.verbose_stream() << "Input: " << input_counter << endl;
		p.verbose_stream() << "Output: " << output_counter << endl;
	}
	return 0;
}
