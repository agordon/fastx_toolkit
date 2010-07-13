#include <istream>
#include <string>
#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "fastx_file.h"
#include "fasta_file.h"
#include "fastq_file.h"

using namespace std;

char peek_first_char(input_stream_wrapper &w)
{
	istream &is(w.stream());
	char c = is.get();
	if (!is) {
		if (is.eof()) {
			cerr << "Input error: File is empty" << endl;
		}
		else {
			cerr << "Input error: failed to read first character from '" << w.filename()
				<< "': " << string_error(errno) << endl;
		}
		exit(1);
	}
	is.unget();

	return c;
}

ISequenceReader* create_fastx_reader(const std::string& filename, int ASCII_quality_offset)
{
	input_stream_wrapper w(filename);

	char c = peek_first_char(w);
	switch(c)
	{
	case '@':
		return new FastqFileReader(w, ASCII_quality_offset);
	case '>':
		return new FastaFileReader(w);

	default:
		cerr << "Error: unknown type (FASTA or FASTQ?) in file '" << filename
			<< "', first character = '" << c << "' (ASCII " << (int)c << ")" << endl;
		exit(1);
	}
}

ISequenceReaderPE* create_fastx_pe_reader(const std::string& filename1, const std::string& filename2, int ASCII_quality_offset)
{
	input_stream_wrapper w1(filename1);
	char c1 = peek_first_char(w1);

	input_stream_wrapper w2(filename2);
	char c2 = peek_first_char(w2);

	if ( c1=='@' && c2=='@' ) {
		return new PE_FastqFileReader(w1,w2,ASCII_quality_offset);
	} else
	if ( c1=='>' && c2=='>' ) {
		return new PE_FastaFileReader(w1,w2);
	}
	//do some extra error checking, trying to find out what's going on
	//give the user a detailed error message
	else if ( c1=='@' && c2=='>' ) {
		cerr << "Input error: input file '" << filename1
			<< "' is a FASTQ file, but input file '" << filename2
			<< "' is a FASTA file. Both input files must be of the same type."
			<< endl;
		exit(1);
	}
	else if ( c1=='>' && c2=='@' ) {
		cerr << "Input error: input file '" << filename1
			<< "' is a FASTA file, but input file '" << filename2
			<< "' is a FASTQ file. Both input files must be of the same type."
			<< endl;
		exit(1);
	}
	else if ( c1!='@' && c1 !='>' ) {
		cerr << "Input error: input file '" << filename1
			<< "' is not a valid FASTA or FASTQ file" << endl;
		exit(1);
	}
	else if ( c2!='@' && c2 !='>' ) {
		cerr << "Input error: input file '" << filename2
			<< "' is not a valid FASTA or FASTQ file" << endl;
		exit(1);
	} else {
		cerr << "Internal error: this should never happen ("
			<< __FILE__ << ":" << __LINE__ << ")" << endl;
		exit(1);
	}
}
