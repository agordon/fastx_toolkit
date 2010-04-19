#ifndef __KMER_CODING__
#define __KMER_CODING__

#include <algorithm>
#include <string>

typedef size_t kmer_type;

extern char kmer_bits_to_char[8];
extern kmer_type kmer_char_to_bit[256];

//Make sure the given k-mer size fits in a "size_t" variable
//Aborts with an error message if not.
//For 32-bit systems,  max size = 32/3 = 10
//For 64-bit systeme,  max size = 64/3 = 21
void assert_valid_kmer_size(size_t size);

inline
kmer_type encode_kmer_string(const std::string& kmer_string)
{
	kmer_type rc=0;
	for (size_t i=0;i<kmer_string.length();++i) {
		rc <<= 3;
		rc |= kmer_char_to_bit[(size_t)kmer_string[i]];
	}
	return rc;
}

inline
std::string decode_kmer_bits(kmer_type kmer_bits)
{
	std::string rc;
	while (kmer_bits>0) {
		kmer_type v = kmer_bits & 0x7 ;
		rc += kmer_bits_to_char[v];
		kmer_bits >>= 3;
	}
	std::reverse(rc.begin(), rc.end());
	return rc;
}

//Runs a simple validation test
//on the coding/decoding routines
void assert_kmer_coding_decoding();

std::string debug_kmer_bitstring(kmer_type bits);
std::string debug_kmer_bits(kmer_type bits);

#endif
