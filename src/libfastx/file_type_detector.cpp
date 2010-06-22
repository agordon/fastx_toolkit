#include <err.h>
#include <ctype.h>
#include <string>
#include <fstream>

#include "file_type_detector.h"

#include <gtextutils/generic_input_stream.h>

using namespace std;

bool file_type_is_fasta_single_line(const std::string &filename)
{
	vector<string> lines;
	file_type_peek_lines(filename, 3, lines);

	//have enough lines ?
	if (lines.size()<2)
		return false;

	if (!is_fasta_id_string(lines[0]))
		return false;
	if (!is_nucleotide_string(lines[1]))
		return false;

	//If the file contains3 or more lines, the third line must be a new
	//sequence id (otherwise, it might be a multi-line FASTA file
	if (lines.size()==3 && !is_fasta_id_string(lines[2]))
		return false;

	return true;
}

bool file_type_is_fastq(const std::string &filename)
{
	vector<string> lines;
	file_type_peek_lines(filename, 5, lines);

	//have enough lines ?
	if (lines.size()<4)
		return false;

	if (!is_fastq_id1_string(lines[0]))
		return false;
	if (!is_nucleotide_string(lines[1]))
		return false;
	if (!is_fastq_id2_string(lines[2]))
		return false;
	//Fixme:  verify quality-score line
	//        (but must allow ASCII and numeric quality values)

	//If the file contains 5 or more lines, the fifth line must be a new sequence id
	if (lines.size()==5 && !is_fastq_id1_string(lines[4]))
		return false;

	return true;
}

bool file_type_is_fastx(const std::string &filename)
{
	return file_type_is_fasta_single_line(filename) || file_type_is_fastq(filename);
}

bool file_type_is_SAM(const std::string &)
{
	errx(1,"SAM Support not implemented yet");
	return false;
}

bool file_type_is_BAM(const std::string &)
{
	errx(1,"BAM Support not implemented yet");
	return false;
}

void file_type_peek_lines(const std::string &filename, size_t num_lines, std::vector< std::string >& /*output*/ lines)
{
	size_t line_count=0;

	generic_input_stream input(filename) ;

	lines.clear();
	lines.reserve(num_lines);

	while ( input && line_count < num_lines ) {
		std::string line;
		getline ( input, line ) ;
		lines.push_back(line);
		++line_count;;
	}
}

bool is_printable_string(const std::string &line)
{
	for (size_t i=0;i<line.size();++i) {
		const char c = line[i];
		if ( ! (isprint(c) || isblank(c)) )
			return false;
	}
	return true;
}

bool is_nucleotide_string(const std::string &line)
{
	if (line.empty())
		return false;
	for (size_t i=0;i<line.size();++i) {
		const char c = toupper(line[i]);
		if ( ! ((c=='A') || (c=='C') || (c=='G') || (c=='T') || (c=='N')) )
			return false;
	}
	return true;
}

bool is_fasta_id_string(const std::string &line)
{
	if (line.size()<2)
		return false;
	if (line.at(0) != '>')
		return false;
	return is_printable_string(line);
}

bool is_fastq_id1_string(const std::string &line)
{
	if (line.empty())
		return false;
	if (line.at(0) != '@')
		return false;
	return is_printable_string(line);
}

bool is_fastq_id2_string(const std::string &line)
{
	if (line.empty())
		return false;
	if (line.at(0) != '+')
		return false;
	if (line.size()==1) //on some machines, the second ID in a FASTQ file can be empty
		return true;
	return is_printable_string(line);
}

bool file_type_is_readable(const std::string &filename)
{
	ifstream f;
	f.open(filename.c_str(),ifstream::in);
	return !f.fail();
}

