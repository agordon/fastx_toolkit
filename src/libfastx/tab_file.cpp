#include <string>
#include <iostream>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>
#include <gtextutils/string_tokenize.h>

#include "sequence.h"
#include "file_type_detector.h"
#include "fastq_file.h"
#include "fasta_file.h"
#include "tab_file.h"

using namespace std;

TabularFileReader::TabularFileReader ( const std::string& filename,int ASCII_Quality_offset  ) :
	_filename(filename), input_stream(filename), line_number(1),
	input_file_format(TAB_FORMAT_UNKNOWN), _ASCII_quality_offset(ASCII_Quality_offset),
	have_cached_line(false)
{
	detect_file_format();
}

TabularFileReader::TabularFileReader ( input_stream_wrapper w, int ASCII_Quality_offset ) :
	_filename(w.filename()), input_stream(w), line_number(1),
	input_file_format(TAB_FORMAT_UNKNOWN), _ASCII_quality_offset(ASCII_Quality_offset),
	have_cached_line(false)
{
	detect_file_format();
}

ISequenceWriter* TabularFileReader::create_fastx_writer(const std::string &filename)
{
	switch(input_file_format)
	{
		case TAB_FORMAT_FASTA:
			return new FastaFileWriter(filename);

		case TAB_FORMAT_FASTQ:
			return new FastqFileWriter(filename);

		case TAB_FORMAT_UNKNOWN:
		default:
			cerr << "Internal error: input_file_format=TAB_FORMAT_UNKNOWN in " << __FILE__
				<< ":" << __LINE__ << endl;
			exit(1);

	}
}

ISequenceWriter* TabularFileReader::create_tabular_writer(const std::string &filename)
{
	return new TabularFileWriter(filename, input_file_format);
}

/*
   Read the first line from the file, try to detect file file type
 */
void TabularFileReader::detect_file_format()
{
	have_cached_line = true;

	if (!getline(input_stream, cached_line)) {
		if (input_stream.eof()) {
			cerr << "Input error: Input file is empty '" << _filename
				<< "'" << endl;
		}
		else {
			cerr << "Input error: failed to read line from '" << _filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
		}
		exit(1);
	}

	input_file_format = detect_line_format ( cached_line ) ;

	//We need to detect the input type,
	//based on number and content of column.
	input_file_format = detect_line_format(cached_line);
	if ( input_file_format == TAB_FORMAT_UNKNOWN ) {
		cerr << "Input error: failed to detect file format (FASTA or FASTQ) from first input line in '" << _filename << "'." << endl;
		exit(1);
	}
}

TABULAR_FILE_FORMAT TabularFileReader::detect_line_format(const std::string& line)
{
	string columns[4];

	size_t found_columns = split_tabular_string(line, columns, 4);

	if (found_columns<2)
		return TAB_FORMAT_UNKNOWN;

	if (found_columns==2) {
		//Assume it's FASTA, just verify
		if (is_printable_string(columns[0]) &&
			is_nucleotide_string(columns[1]) )
			return TAB_FORMAT_FASTA;

		return TAB_FORMAT_UNKNOWN;
	}

	if (found_columns==4) {
		//Assume it's FASTQ, just verify
		if (is_printable_string(columns[0]) &&
			is_nucleotide_string(columns[1]) &&
			(columns[3].size() == columns[1].size()) &&
			!is_nucleotide_string(columns[3]))
			return TAB_FORMAT_FASTQ;

		return TAB_FORMAT_UNKNOWN;
	}
	return TAB_FORMAT_UNKNOWN;
}

bool TabularFileReader::read_next_sequence(Sequence& output)
{
	string line;

	output.clear();

	if (have_cached_line) {
		line = cached_line ;
		have_cached_line = false;
	}
	else {
		if (!getline(input_stream, line)) {
			if (input_stream.eof())
				return false;

			cerr << "Input error: failed to read line from '" << _filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}
	}

	string columns[4];
	split_tabular_string(line, columns, 4);

	//Verify input
	if (!is_printable_string(columns[0])) {
			cerr << "Input error: Invalid Sequence-ID (column 1 = \"" << columns[0] << "\") in '" << _filename
				<< "' line " << line_number << endl;
			exit(1);
	}
	if (!is_nucleotide_string(columns[1])) {
			cerr << "Input error: Invalid Nucleotide value (column 2 = \"" << columns[1] << "\") in '" << _filename
				<< "' line " << line_number << endl;
			exit(1);
	}

	if ( input_file_format == TAB_FORMAT_FASTQ ) {
		if (columns[3].size() != columns[1].size()) {
				cerr << "Input error: Invalid Quality-score value (column 4 = \"" << columns[3] << "\", does not match length of column 1) in '" << _filename << "' line " << line_number << endl;
				exit(1);
		}
	}

	// Set values
	output.id = columns[0];
	output.nucleotides = columns[1] ;
	if ( input_file_format == TAB_FORMAT_FASTQ ) {
		output.id2 = columns[2] ;
		FastqFileReader::convert_ascii_quality_score_line ( columns[3], output.quality, _ASCII_quality_offset ) ;
		output.ASCII_quality_offset = _ASCII_quality_offset ;
		output.ASCII_quality_scores = true;
	}

	++line_number;

	return true;
}

TabularFileWriter::TabularFileWriter ( const std::string& filename, TABULAR_FILE_FORMAT _output_format ) :
	_filename ( filename.empty()?"stdout":filename ) ,
	output_stream ( filename ), output_file_format(_output_format)
{
}

void TabularFileWriter::write_sequence(const Sequence& seq)
{
	if (seq.id.empty()) {
		cerr << "Internal error: about to write an empty sequence ID." << endl;
		exit(1);
	}

	output_stream << seq.id << "\t" << seq.nucleotides;
	if ( output_file_format == TAB_FORMAT_FASTA ) {
	} else
	if ( output_file_format == TAB_FORMAT_FASTQ ) {
		output_stream << "\t" << seq.id2 << "\t" ;

		for ( size_t i=0;i<seq.quality.size();++i) {
			const char c = (char)(seq.quality[i] + seq.ASCII_quality_offset);
			output_stream << c;
		}
	} else {
		cerr << "Internal error: output file format is unknown (FASTA or FASTQ?) in "
			<< __FILE__ << ":" << __LINE__ << endl;
		exit(1);
	}

	output_stream << endl;

	if (!output_stream) {
		cerr << "Output error: failed to write data to '" << _filename
			<< "': " << string_error(errno) << endl;
		exit(1);
	}
}

/*
   Paired-end reader
 */

PE_TabularFileReader::PE_TabularFileReader ( const std::string& filename, int ASCII_quality_offset ) :
	_filename(filename), input_stream(filename), line_number(1),
	input_file_format(TAB_FORMAT_UNKNOWN), _ASCII_quality_offset(ASCII_quality_offset),
	have_cached_line(false)
{
	detect_file_format();
}

PE_TabularFileReader::PE_TabularFileReader ( input_stream_wrapper w, int ASCII_quality_offset ) :
	_filename(w.filename()), input_stream(w), line_number(1),
	input_file_format(TAB_FORMAT_UNKNOWN), _ASCII_quality_offset(ASCII_quality_offset),
	have_cached_line(false)
{
	detect_file_format();
}

void PE_TabularFileReader::detect_file_format()
{
	have_cached_line = true;

	if (!getline(input_stream, cached_line)) {
		if (input_stream.eof()) {
			cerr << "Input error: Input file is empty '" << _filename
				<< "'" << endl;
		}
		else {
			cerr << "Input error: failed to read line from '" << _filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
		}
		exit(1);
	}

	input_file_format = detect_line_format ( cached_line ) ;

	//We need to detect the input type,
	//based on number and content of column.
	input_file_format = detect_line_format(cached_line);
	if ( input_file_format == TAB_FORMAT_UNKNOWN ) {
		cerr << "Input error: failed to detect file format (PE-FASTA or PE-FASTQ) from first input line in '" << _filename << "'." << endl;
		exit(1);
	}
}

bool PE_TabularFileReader::valid_pe_fasta_line(const string columns[4])
{
	return	(is_printable_string(columns[0]) &&
			is_nucleotide_string(columns[1]) &&
			is_printable_string(columns[2]) &&
			is_nucleotide_string(columns[3]));
}

bool PE_TabularFileReader::valid_pe_fastq_line(const string columns[8])
{
	return (is_printable_string(columns[0]) &&
		    is_nucleotide_string(columns[1]) &&
		    (columns[3].size() == columns[1].size()) &&
		    (!is_nucleotide_string(columns[3]))

			&&
		     is_printable_string(columns[4]) &&
		    is_nucleotide_string(columns[5]) &&
		    (columns[7].size() == columns[5].size()) &&
		    (!is_nucleotide_string(columns[7]))
		   );
}

TABULAR_FILE_FORMAT PE_TabularFileReader::detect_line_format(const std::string& line)
{
	string columns[4];

	size_t found_columns = split_tabular_string(line, columns, 4);

	if (found_columns<4)
		return TAB_FORMAT_UNKNOWN;

	if (found_columns==4 && valid_pe_fasta_line(columns))
		return TAB_FORMAT_FASTA;

	if (found_columns==8 && valid_pe_fastq_line(columns))
		return TAB_FORMAT_FASTQ;

	return TAB_FORMAT_UNKNOWN;
}

void PE_TabularFileReader::parse_pe_fasta_line(const std::string& line, Sequence& out_seq1, Sequence& out_seq2)
{
	string columns[4];
	int num_columns = split_tabular_string(line, columns, 4);

	if (num_columns!=4) {
		cerr << "Input error: expecting 4 fields for PE-FASTA tabular file in '"
			<< _filename << "' line " << line_number
		       << ", found " << num_columns << " fields ("
		      << line << ")" << endl;
		exit(1);
	}
	if (!valid_pe_fasta_line(columns)) {
		cerr << "Input error: Invalid values for PE-FASTA tabular file in '"
			<< _filename << "' line " << line_number
			<< " (" << line << ")" << endl;
		exit(1);
	}

	out_seq1.id = columns[0];
	out_seq1.nucleotides = columns[1];
	out_seq1.ASCII_quality_offset = _ASCII_quality_offset ;
	out_seq1.ASCII_quality_scores = true;

	out_seq2.id = columns[2];
	out_seq2.nucleotides = columns[3];
	out_seq2.ASCII_quality_offset = _ASCII_quality_offset ;
	out_seq2.ASCII_quality_scores = true;
}

void PE_TabularFileReader::parse_pe_fastq_line(const std::string& line, Sequence& out_seq1, Sequence& out_seq2)
{
	string columns[8];
	int num_columns = split_tabular_string(line, columns, 8);

	if (num_columns!=8) {
		cerr << "Input error: expecting 8 fields for PE-FASTA tabular file in '"
			<< _filename << "' line " << line_number
		       << ", found " << num_columns << " fields ("
		      << line << ")" << endl;
		exit(1);
	}
	if (!valid_pe_fastq_line(columns)) {
		cerr << "Input error: Invalid values for PE-FASTQ tabular file in '"
			<< _filename << "' line " << line_number
			<< " (" << line << ")" << endl;
		exit(1);
	}

	out_seq1.id = columns[0];
	out_seq1.nucleotides = columns[1];
	out_seq1.id2 = columns[2];
	FastqFileReader::convert_ascii_quality_score_line ( columns[3], out_seq1.quality, _ASCII_quality_offset ) ;
	out_seq1.ASCII_quality_offset = _ASCII_quality_offset ;
	out_seq1.ASCII_quality_scores = true;

	out_seq2.id = columns[4];
	out_seq2.nucleotides = columns[5];
	out_seq2.id2 = columns[6];
	FastqFileReader::convert_ascii_quality_score_line ( columns[7], out_seq2.quality, _ASCII_quality_offset ) ;
	out_seq2.ASCII_quality_offset = _ASCII_quality_offset ;
	out_seq2.ASCII_quality_scores = true;
}

bool PE_TabularFileReader::read_next_sequence(Sequence& out_seq1, Sequence& out_seq2)
{
	string line;

	out_seq1.clear();
	out_seq2.clear();

	if (have_cached_line) {
		line = cached_line ;
		have_cached_line = false;
	}
	else {
		if (!getline(input_stream, line)) {
			if (input_stream.eof())
				return false;

			cerr << "Input error: failed to read line from '" << _filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}
	}


	switch(input_file_format)
	{
	case TAB_FORMAT_FASTA:
		parse_pe_fasta_line(line, out_seq1, out_seq2);
		break;

	case TAB_FORMAT_FASTQ:
		parse_pe_fastq_line(line,out_seq2, out_seq2);
		break;

	default:
	case TAB_FORMAT_UNKNOWN:
		cerr << "Internal error: input_file_format=TAB_FORMAT_UNKNOWN in "
			<< __FILE__ << ":" << __LINE__ << endl;
		exit(1);
		break;
	}



	++line_number;

	return true;
}

ISequenceWriterPE* PE_TabularFileReader::create_fastx_writer(const std::string& filename1, const std::string &filename2)
{
	switch(input_file_format)
	{
		case TAB_FORMAT_FASTA:
			return new PE_FastaFileWriter(filename1, filename2);

		case TAB_FORMAT_FASTQ:
			return new PE_FastqFileWriter(filename1, filename2);

		case TAB_FORMAT_UNKNOWN:
		default:
			cerr << "Internal error: input_file_format=TAB_FORMAT_UNKNOWN in " << __FILE__
				<< ":" << __LINE__ << endl;
			exit(1);

	}
}

ISequenceWriterPE* PE_TabularFileReader::create_tabular_writer(const std::string& filename)
{
	return new PE_TabularFileWriter(filename, input_file_format);
}


/*
   Paired-end writer
 */

PE_TabularFileWriter::PE_TabularFileWriter ( const std::string& filename , TABULAR_FILE_FORMAT _output_format ) :
	_filename(filename), output_stream(filename), output_file_format(_output_format)
{
}

void PE_TabularFileWriter::write_tabular_sequence(std::ostream& ostm, const Sequence& seq)
{
	ostm << seq.id << "\t" << seq.nucleotides;

	switch(output_file_format)
	{
	case TAB_FORMAT_FASTA:
		break;
	case TAB_FORMAT_FASTQ:
		output_stream << "\t" << seq.id2 << "\t" ;
		for ( size_t i=0;i<seq.quality.size();++i) {
			const char c = (char)(seq.quality[i] + seq.ASCII_quality_offset);
			output_stream << c;
		}
		break;

	default:
	case TAB_FORMAT_UNKNOWN:
		cerr << "Internal error: output file format is unknown (FASTA or FASTQ?) in "
			<< __FILE__ << ":" << __LINE__ << endl;
		exit(1);
		break;
	}
}

void PE_TabularFileWriter::write_sequence(const Sequence& seq1, const Sequence& seq2)
{
	if (seq1.id.empty() || seq2.id.empty()) {
		cerr << "Internal error: about to write an empty sequence ID." << endl;
		exit(1);
	}

	write_tabular_sequence(output_stream, seq1);
	output_stream << "\t";
	write_tabular_sequence(output_stream, seq2);
	output_stream << endl;

	if (!output_stream) {
		cerr << "Output error: failed to write data to '" << _filename
			<< "': " << string_error(errno) << endl;
		exit(1);
	}
}
