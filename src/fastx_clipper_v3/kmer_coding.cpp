#include <assert.h>
#include <ctype.h>
#include <err.h>
#include <iostream>
#include <sstream>
#include "kmer_coding.h"

char kmer_bits_to_char[8] =
{
	'*', //0 - not a valid char
	'A', //1 - A
	'C', //2 - C
	'G', //3 - G
	'T', //4 - T
	'N', //5 - N
	'*', //6 - not a valid char
	'*', //7 - not a valid char
};

size_t kmer_char_to_bit[256] =
{
	0,	0,	0,	0,	0,	0,	0,	0, //0
	0,	0,	0,	0,	0,	0,	0,	0, //8
	0,	0,	0,	0,	0,	0,	0,	0, //16
	0,	0,	0,	0,	0,	0,	0,	0, //24
	0,	0,	0,	0,	0,	0,	0,	0, //32
	0,	0,	0,	0,	0,	0,	0,	0, //40
	0,	0,	0,	0,	0,	0,	0,	0, //48
	0,	0,	0,	0,	0,	0,	0,	0, //56
	0,	1,	0,	2,	0,	0,	0,	3, //64
	0,	0,	0,	0,	0,	0,	5,	0, //72
	0,	0,	0,	0,	4,	0,	0,	0, //80
	0,	0,	0,	0,	0,	0,	0,	0, //88
	0,	1,	0,	2,	0,	0,	0,	3, //96
	0,	0,	0,	0,	0,	0,	5,	0, //104
	0,	0,	0,	0,	4,	0,	0,	0, //112
	0,	0,	0,	0,	0,	0,	0,	0, //120
	0,	0,	0,	0,	0,	0,	0,	0, //128
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0
};

using namespace std;

string debug_kmer_bitstring(kmer_type bits)
{
	stringstream ss;

	while (bits>0) {
		kmer_type v = bits & 0x7 ;

		ss << v << "/" << kmer_bits_to_char[ v ] << " " ;
		bits >>= 3;
	}

	return ss.str();
}

string debug_kmer_bits(kmer_type bits)
{
	stringstream ss;

	while (bits>0) {
		kmer_type v = bits & 0x1 ;

		ss << v ;
		bits >>= 1;
	}

	return ss.str();
}

void recursive_check_kmer_decoding(std::string &input, unsigned int pos)
{
	if (pos==0) {
		//next permutation done - check the coding/decoding
		kmer_type v = encode_kmer_string(input);
		std::string decoded_value = decode_kmer_bits(v);

		cerr << "Input = " << input
			<< "\tValue = " << hex << v << oct
			<< "\tBits = " << debug_kmer_bits(v)
			<< "\tBits2 = " << debug_kmer_bitstring(v)
			<< "\tChar = " << decoded_value
			<< endl;

		assert(decoded_value == input);
		return ;
	}

	static const char *nucleotide = "ACGTN" ;
	for (size_t i=0;i<5;++i) {
		input[pos-1] = nucleotide[i];
		recursive_check_kmer_decoding(input, pos-1);
	}
}

void assert_kmer_coding_decoding()
{
	std::string s("ACGTNacgtn");

	//Validate translation tables
	for (size_t i=0;i<s.length();++i) {
		kmer_type v = kmer_char_to_bit[ (size_t) s[i] ];
		assert(v!=0);

		char c = kmer_bits_to_char[v];

		/* cerr << "input char = " << s[i]
			<< " kmer_char_to_bits = " << std::hex << v << std::oct
			<< " kmet_bits_to_char = " << c<< std::endl;
		*/
		assert(toupper(s[i]) == c);
	}


	//Validate encoding/decoding functions
	string data = "12345678";
	recursive_check_kmer_decoding(data,data.length());
}

void assert_valid_kmer_size(size_t size)
{
	if ((sizeof(size_t)*8/3) >= size)
		return;
	errx(1,"K-mer size of %zu can not be used on this system (due to bit-size of 'size_t' = %zu bits). Maximum k-mer size is %zu",
			size, sizeof(size_t)*8, (sizeof(size_t)*8)/3);
}
