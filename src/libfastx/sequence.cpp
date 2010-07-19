#include <iostream>
#include <cstdlib>
#include "sequence.h"

using namespace std;

Sequence::Sequence()
{
	clear();
}

void Sequence::clear()
{
	id.clear();
	nucleotides.clear();
	id2.clear();
//	quality.clear();
	quality_cached_line.clear();

	ASCII_quality_offset = 0;
	ASCII_quality_scores = true;
}

void Sequence::convert_ascii_quality_score_line ( const std::string&  )
{
	/*
	quality.clear();
	quality.resize(quality_line.length());
	for ( size_t i = 0 ; i< quality_line.length(); ++i ) {
		int value = ((int)quality_line[i]) - ASCII_quality_offset;
		quality[i] = value;
	}
	*/
}

void Sequence::convert_numeric_quality_score_line ( const std::string &numeric_quality_line, const std::string& _filename, size_t line_number)
{
	size_t index;
	const char *quality_tok;
	char *endptr;
	int quality_value;

	std::vector<int> quality;

	quality.clear();
	quality.resize(nucleotides.length());

	index=0;
	quality_tok = numeric_quality_line.c_str();
	do {
		//read the quality score as an integer value
		quality_value = strtol(quality_tok, &endptr, 10);
		if (endptr == quality_tok) {
			cerr << "Input error: invalid quality score data (" << quality_tok << ") in '" << _filename << "' line " << line_number << endl ;
			exit(1);
		}
		if (index>=quality.size()) {
			cerr << "Input error: too many quality score values in '" << _filename << "' line " << line_number
			<< ": expected " << nucleotides.size() << " numbers."	<< endl ;
			exit(1);
		}

		//convert it ASCII (as per solexa's encoding)
		quality[index] = ( quality_value - ASCII_quality_offset) ;

		index++;
		quality_tok = endptr;
	} while (quality_tok != NULL && *quality_tok!='\0') ;
}

