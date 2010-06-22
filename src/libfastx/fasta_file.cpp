#include <string>
#include <iostream>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "file_type_detector.h"
#include "fasta_file.h"

using namespace std;

FastaFileReader::FastaFileReader ( const std::string& filename ) :
	_filename(filename), input_stream(filename), line_number(0)
{
}

bool FastaFileReader::read_next_sequence(Sequence& output)
{
	string id;

	output.clear();

	if (!getline(input_stream, id))
		return false;	//EOF

	if (!is_fasta_id_string(id)) {
		cerr << "Input error: Invalid FASTA ID value (" << id << ") in '" << _filename
			<< "' line " << line_number << endl;
			exit(1);
	}
	++line_number;

	const string &id_no_prefix ( id.substr(1) ) ;
	output.id = id_no_prefix;

	string nuc;
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

	output.nucleotides = nuc ;

	return true;
}

FastaFileWriter::FastaFileWriter ( const std::string& filename ) :
	_filename ( filename ) ,
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
}
