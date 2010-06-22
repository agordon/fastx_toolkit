#include <string>
#include <iostream>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "file_type_detector.h"
#include "fasta_file.h"

using namespace std;

FastaFileReader::FastaFileReader ( const std::string& filename ) :
	_filename(filename), input_stream(filename), line_number(1)
{
}

FastaFileReader::FastaFileReader ( input_stream_wrapper w ) :
	_filename(w.filename()), input_stream(w), line_number(1)
{
}

ISequenceWriter* FastaFileReader::create_writer(const std::string &filename)
{
	return new FastaFileWriter(filename);
}

bool FastaFileReader::read_next_sequence(Sequence& output)
{
	string id;

	output.clear();

	if (!getline(input_stream, id)) {
		if (input_stream.eof())
			return false;

		cerr << "Input error: failed to read ID line from '" << _filename
			<< "' line " << line_number << ":" << string_error(errno) << endl;
		exit(1);
	}

	if (!is_fasta_id_string(id)) {
		if (is_nucleotide_string(id)) {
			cerr << "Input error: Found multi-line nucleotides sequence in '" << _filename
				<< "' line " << line_number << ". This program can only read single-line FASTA sequences. Use 'fasta_formatter' to convert the file to a single-line FASTA file." << endl;
		}
		else {
			cerr << "Input error: Invalid FASTA ID value (" << id << ") in '" << _filename
				<< "' line " << line_number << endl;
		}
		exit(1);
	}
	++line_number;

	const string &id_no_prefix ( id.substr(1) ) ;
	output.id = id_no_prefix;

	string nuc;
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
		if ((!nuc.empty()) && nuc.at(0)=='>') {
			cerr << "Input error: Empty sequence found in '" << _filename
				<< "' line " << (line_number-1) << " (two consecutive >ID lines) " << endl;
		}
		else {
			cerr << "Input error: Invalid nucleotides line (" << nuc << ") in '" << _filename
				<< "' line " << line_number << endl;
		}
		exit(1);
	}
	++line_number;

	output.nucleotides = nuc ;

	return true;
}

FastaFileWriter::FastaFileWriter ( const std::string& filename ) :
	_filename ( filename.empty()?"stdout":filename ) ,
	output_stream ( filename )
{
}

void FastaFileWriter::write_sequence(const Sequence& seq)
{
	if (seq.id.empty()) {
		cerr << "Internal error: about to write an empty sequence ID." << endl;
		exit(1);
	}

	output_stream << ">" << seq.id << endl;
	output_stream << seq.nucleotides << endl;
	if (!output_stream) {
		cerr << "Output error: failed to write data to '" << _filename
			<< "': " << string_error(errno) << endl;
		exit(1);
	}
}
