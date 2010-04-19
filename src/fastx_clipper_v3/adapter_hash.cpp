#include <err.h>
#include <iostream>
#include <string>
#include <tr1/unordered_set>

#include "adapter_hash.h"

using namespace std;

SequenceHash::SequenceHash()
	: kmer_size(0)
{
}

void SequenceHash::set_sequence(const std::string& _sequence, size_t _kmer_size)
{
	sequence = _sequence ;
	kmer_size = _kmer_size;

	build_hash();
}

void SequenceHash::build_hash()
{
	int k_step = 1;
	tr1::unordered_set<string> kmers_with_multiple_offsets;

	for (size_t  i=0;i<sequence.length() - kmer_size ; i+=k_step) {

		const string kmer = sequence.substr(i, kmer_size);

		kmers_hash::iterator it = kmers.find(kmer);
		if ( it == kmers.end() ) {
			if ( kmers_with_multiple_offsets.find(kmer)==kmers_with_multiple_offsets.end()) {
				//This is a new kmer for this sequence - add it with one offset
				kmers.insert ( kmers_hash::value_type(kmer, i) );
			} else {
				//This k-Mer is marked as having multiple offsets - don't add it again.
			}
		}
		else {
			//This k-mer was already found in the hash - it has more than one offset.
			//Mark it and delete it
			kmers_with_multiple_offsets.insert ( kmer ) ;
			kmers.erase(it);
		}
	}


	//Warning to the user:
	//if the first k-mer has multiple offsets (and was removed)
	//it will make it hard to detect some extreme cases of the adapter
	string first_kmer = sequence.substr(0,kmer_size);
	if (kmers.find(first_kmer)==kmers.end()) {
		cerr << "Warning: the first k-mer (" << first_kmer << ") in the adapter sequence appears multiple times (in the adapter)." << endl
		       << "   This will make adapter detection worse. Some cases might be mis-detected. " << endl
		       << "   Consider using a trimmed adapter, skipping the first few nucleotdes," << endl
		       << "   or a longer k-mer size (-k)" << endl;
	}
}

void SequenceHash::debug_print_hash(std::ostream &strm)
{
	kmers_hash::iterator it = kmers.begin();
	while ( it != kmers.end() ) {
		const string &kmer_sequence = it->first;
		const size_t &data = it->second;

		strm << "kmer = " << kmer_sequence ;
		strm << " offset = " << data;
		strm << endl;
		++it;
	}
}
