#ifndef __ADAPTER_HASH__
#define __ADAPTER_HASH__

#include <vector>
#include <tr1/unordered_map>
#include <string>
#include <limits.h>
#include <ostream>
#include <stdlib.h>

typedef std::string kmer_type;
typedef std::tr1::unordered_map<kmer_type, size_t> kmers_hash;

class SequenceHash
{
private:
	std::string sequence;
	size_t kmer_size;

	kmers_hash kmers;

public:
	SequenceHash();

	void set_sequence(const std::string& _sequence, size_t _kmer_size);
	void debug_print_hash(std::ostream &strm);

	size_t get_kmer_size() const { return kmer_size; }

	bool kmer_exists(const kmer_type& kmer) const
	{
		return kmers.find(kmer) != kmers.end();
	}

	size_t get_kmer_offset(const kmer_type& kmer) const
	{
		return kmers.find(kmer)->second;
	}

private:
	void build_hash();
};

#endif

