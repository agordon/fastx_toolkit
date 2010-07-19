#include <string>
#include <iostream>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "file_type_detector.h"
#include "fastq_file.h"
#include "tab_file.h"

using namespace std;

FastqFileReader::FastqFileReader ( const std::string& filename, int ASCII_quality_offset ) :
	_filename(filename), input_stream(filename), line_number(1), _ASCII_quality_offset(ASCII_quality_offset)
{
}

FastqFileReader::FastqFileReader ( input_stream_wrapper w, int ASCII_quality_offset ):
	_filename(w.filename()), input_stream(w), line_number(1), _ASCII_quality_offset(ASCII_quality_offset)
{
}

ISequenceWriter* FastqFileReader::create_fastx_writer(const std::string &filename)
{
	return new FastqFileWriter(filename);
}

ISequenceWriter* FastqFileReader::create_tabular_writer(const std::string &filename)
{
	return new TabularFileWriter(filename, TAB_FORMAT_FASTQ);
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

		cerr << "Input error: failed to read ID line from '" << _filename
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
		if (input_stream.eof()) {
			cerr << "Input error: premature End-of-File in '" << _filename
				<< "' line " << line_number << ": expecting nucleotides line." << endl;
		}
		else {
			cerr << "Input error: failed to read nucleotides from '" << _filename
				<< "' line " << line_number << ": " << string_error(errno) << endl;
		}
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
		if (input_stream.eof()) {
			cerr << "Input error: premature End-of-File in '" << _filename
				<< "' line " << line_number << ": expecting ID line." << endl;
		}
		else {
			cerr << "Input error: failed to read ID-2 from '" << _filename
				<< "' line " << line_number << ": " << string_error(errno) << endl;
		}
		exit(1);
	}
	if (!is_fastq_id2_string(id2)) {
		cerr << "Input error: Invalid FASTQ ID-2 value (" << id << ") in '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	const string &id2_no_prefix ( id2.substr(1) ) ;
	++line_number;
	output.id2 = id2_no_prefix;

	//Line 4 - quality scores
	if (!getline(input_stream, quality)) {
		if (input_stream.eof()) {
			cerr << "Input error: premature End-of-File in '" << _filename
				<< "' line " << line_number << ": expecting quality-score line." << endl;
		}
		else {
			cerr << "Input error: failed to read quality-scores from '" << _filename
				<< "' line " << line_number << ": " << string_error(errno) << endl;
		}
		exit(1);
	}

	/*
	if ( quality.length() == nuc.length() ) {
		//ASCII quality scores
		convert_ascii_quality_score_line ( quality, output.quality, _ASCII_quality_offset ) ;
		output.ASCII_quality_scores = true ;
	} else {
		//Numeric Quality Scores
		convert_numeric_quality_score_line ( quality, output.quality ) ;
		output.ASCII_quality_scores = false ;
	}
	*/
	output.quality_cached_line = quality ;
	output.ASCII_quality_offset = _ASCII_quality_offset ;

	/*
	if (output.quality.size() != output.nucleotides.size()) {
		cerr << "Input error: number of quality-score values (" << output.quality.size()
			<< ") doesn't match number of nuecltoes (" << output.nucleotides.size()
			<< ") in file '" << _filename
			<< "' line " << line_number << endl;
		exit(1);
	}
	*/
	++line_number;

	return true;
}


FastqFileWriter::FastqFileWriter ( const std::string& filename ) :
	_filename ( filename.empty()?"stdout":filename ) ,
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
	if (seq.quality_cached_line.empty()) {
		if (seq.id2.empty())
			output_stream << "+" << seq.id << endl;
		else
			output_stream << "+" << seq.id2 << endl;

		for (size_t i=0;i<seq.nucleotides.length();++i)
			output_stream << "b";
		output_stream << endl;

	} else {
		output_stream << "+" << seq.id2 << endl;
		output_stream << seq.quality_cached_line << endl ;
	}
	if (!output_stream) {
		cerr << "Output error: failed to write data to '" << _filename
			<< "': " << string_error(errno) << endl;
		exit(1);
	}
}


PE_FastqFileReader::PE_FastqFileReader ( const std::string& filename1, const std::string& filename2,
	       int ASCII_quality_offset ) :
	end1(filename1, ASCII_quality_offset), end2(filename2, ASCII_quality_offset)
{
}

PE_FastqFileReader::PE_FastqFileReader ( input_stream_wrapper file1, input_stream_wrapper file2, int ASCII_quality_offset ) :
	end1(file1, ASCII_quality_offset), end2(file2, ASCII_quality_offset)
{
}

bool PE_FastqFileReader::read_next_sequence(Sequence& /*output*/ seq1, Sequence& /*output*/ seq2)
{
	bool b1 = end1.read_next_sequence(seq1);
	bool b2 = end2.read_next_sequence(seq2);

	if (b1 != b2) {
		cerr << "Input error: Paired-end FASTQ file mismatch: file '" <<
			(b1 ? end2.filename() : end1.filename()) << "' ended before file '" <<
			(b1 ? end1.filename() : end2.filename()) << "'. This program requires both FASTA file to have the same number of sequences." << endl;
		exit(1);
	}
	return b1;
}

ISequenceWriterPE* PE_FastqFileReader::create_fastx_writer(const std::string& filename1, const std::string &filename2)
{
	return new PE_FastqFileWriter(filename1, filename2);
}

ISequenceWriterPE* PE_FastqFileReader::create_tabular_writer(const std::string& filename)
{
	return new PE_TabularFileWriter(filename, TAB_FORMAT_FASTQ);
}

PE_FastqFileWriter::PE_FastqFileWriter(const std::string& filename1, const std::string& filename2) :
	end1(filename1), end2(filename2)
{
}

void PE_FastqFileWriter::write_sequence(const Sequence& seq1, const Sequence& seq2)
{
	end1.write_sequence(seq1);
	end2.write_sequence(seq2);
}

