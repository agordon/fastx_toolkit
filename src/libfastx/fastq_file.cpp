#include <string>
#include <iostream>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "file_type_detector.h"
#include "fastq_file.h"

using namespace std;

FastqFileReader::FastqFileReader ( const std::string& filename, int ASCII_quality_offset ) :
	_filename(filename), input_stream(filename), line_number(0), _ASCII_quality_offset(ASCII_quality_offset)
{
}

bool FastqFileReader::read_next_sequence(Sequence& output)
{
	string id;
	string nuc;
	string id2;
	string quality;

	output.clear();

	//Line 1 - ID
	if (!getline(input_stream, id)) {
		if (input_stream.eof())
			return false;

		cerr << "Input error: failed to read ID from '" << _filename
			<< "' line " << line_number << ":" << string_error(errno) << endl;
		exit(1);
	}

	if (!is_fastq_id1_string(id)) {
		cerr << "Input error: Invalid FASTQ ID value (" << id << ") in '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	const string &id_no_prefix ( id.substr(1) ) ;
	++line_number;
	output.id = id_no_prefix;

	//line 2 - nucleotides
	if (!getline(input_stream, nuc)) {
		cerr << "Input error: failed to read nucleotides from '" << _filename
			<< "' line " << line_number << ":" << string_error(errno) << endl;
		exit(1);
	}

	if (!is_nucleotide_string(nuc)) {
		cerr << "Input error: Invalid nucleotides line (" << nuc << ") in '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	++line_number;
	output.nucleotides = nuc;

	//Line 3 - ID22
	if (!getline(input_stream, id2)) {
		cerr << "Input error: failed to read ID-2 from '" << _filename
			<< "' line " << line_number << ":" << string_error(errno) << endl;
		exit(1);
	}
	if (!is_fastq_id2_string(id2)) {
		cerr << "Input error: Invalid FASTQ ID-2 value (" << id << ") in '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	++line_number;
	output.id2 = id2;

	//Line 4 - quality scores
	if (!getline(input_stream, quality)) {
		cerr << "Input error: failed to read quality-scores from '" << _filename
			<< "' line " << line_number << ":" << string_error(errno) << endl;
		exit(1);
	}

	if ( quality.length() == nuc.length() ) {
		//ASCII quality scores
		convert_ascii_quality_score_line ( quality, output.quality ) ;
		output.ASCII_quality_scores = true ;
	} else {
		//Numeric Quality Scores
		convert_numeric_quality_score_line ( quality, output.quality ) ;
		output.ASCII_quality_scores = false ;
	}
	output.ASCII_quality_offset = _ASCII_quality_offset ;

	if (output.quality.size() != output.nucleotides.size()) {
		cerr << "Input error: number of quality-score values (" << output.quality.size()
			<< ") doesn't match number of nuecltoes (" << output.nucleotides.size()
			<< ") in file '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	++line_number;

	return true;
}

void FastqFileReader::convert_ascii_quality_score_line ( const std::string& quality_line, std::vector<int>& output)
{
	output.clear();
	output.reserve(100);
	for ( size_t i = 0 ; i< quality_line.size(); ++i ) {
		int value = ((int)quality_line[i]) - _ASCII_quality_offset ;
		output.push_back(value);
	}
}

void FastqFileReader::convert_numeric_quality_score_line ( const std::string &numeric_quality_line, std::vector<int>& output)
{
	size_t index;
	const char *quality_tok;
	char *endptr;
	int quality_value;

	output.clear();
	output.reserve(100);

	index=0;
	quality_tok = numeric_quality_line.c_str();
	do {
		//read the quality score as an integer value
		quality_value = strtol(quality_tok, &endptr, 10);
		if (endptr == quality_tok) {
			cerr << "Error: invalid quality score data in file " << _filename << " line " << line_number << endl ;
			exit(1);
		}

		//convert it ASCII (as per solexa's encoding)
		output.push_back( quality_value - _ASCII_quality_offset ) ;

		index++;
		quality_tok = endptr;
	} while (quality_tok != NULL && *quality_tok!='\0') ;
}

FastqFileWriter::FastqFileWriter ( const std::string& filename ) :
	_filename ( filename ) ,
	output_stream ( filename )
{
}

void FastqFileWriter::write_sequence(const Sequence& seq)
{
	if (seq.id.empty()) {
		cerr << "Internal error: about to write an empty sequence ID." << endl;
		exit(1);
	}

	output_stream << "@" << seq.id << endl;
	output_stream << seq.nucleotides << endl;

	//If seq comes from a FASTA file and doesn't have id2/qualities - fake them
	if (seq.quality.empty()) {
		if (seq.id2.empty())
			output_stream << "+" << seq.id << endl;
		else
			output_stream << "+" << seq.id2 << endl;

		if (seq.ASCII_quality_scores) {
			for ( size_t i=0;i<seq.quality.size();++i) {
				const char c = (char)(0 + seq.ASCII_quality_offset);
				output_stream << c;
			}
		} else {
			for ( size_t i=0;i<seq.quality.size();++i) {
				const int val = (0 + seq.ASCII_quality_offset);
				if (i>0)
					output_stream << " " ;
				output_stream << val;
			}
		}
		output_stream << endl;


	} else {
		output_stream << "+" << seq.id2 << endl;

		if (seq.ASCII_quality_scores) {
			for ( size_t i=0;i<seq.quality.size();++i) {
				const char c = (char)(seq.quality[i] + seq.ASCII_quality_offset);
				output_stream << c;
			}
		} else {
			for ( size_t i=0;i<seq.quality.size();++i) {
				const int val = (seq.quality[i] + seq.ASCII_quality_offset);
				if (i>0)
					output_stream << " " ;
				output_stream << val;
			}
		}
		output_stream << endl;
	}
}
