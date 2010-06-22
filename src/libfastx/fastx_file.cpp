#include <istream>
#include <string>
#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "sequence.h"
#include "fastx_file.h"
#include "fasta_file.h"
#include "fastq_file.h"

using namespace std;

ISequenceReader* create_fastx_reader(const std::string& filename, int ASCII_quality_offset)
{
	input_stream_wrapper w(filename);

	istream &is(w.stream());
	char c = is.get();
	is.unget();

	if (c =='@') {
		return new FastqFileReader(w, ASCII_quality_offset);
	} else
	if (c =='>') {
		return new FastaFileReader(w);
	} else {
		cerr << "Error: unknown type (FASTA or FASTQ?) in file '" << filename
			<< "', first character = '" << c << "' (ASCII " << (int)c << ")" << endl;
		exit(1);
	}
}
