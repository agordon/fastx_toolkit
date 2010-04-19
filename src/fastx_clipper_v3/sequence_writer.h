#ifndef __SEQUENCE_WRITER__
#define __SEQUENCE_WRITER__

#include <string>
#include <memory>
#include <gtextutils/stream_wrapper.h>

class SequenceWriter
{
	std::string	filename;
	OutputStreamWrapper stream;
	std::auto_ptr<FastxWriter>	writer;

public:
	enum TYPE {
		FASTQ = 1,
		FASTA
	};
	SequenceWriter(const std::string& _filename, TYPE type, const std::string& suffix) :
		filename (_filename),
		stream ( _filename )
	{
		switch(type)
		{
		case FASTQ:
			writer = std::auto_ptr<FastxWriter> ( new FastqWriter ( stream.stream(), suffix ) ) ;
			break;
		case FASTA:
			writer = std::auto_ptr<FastxWriter> ( new FastaWriter ( stream.stream(), suffix ) ) ;
			break;
		default:
			errx(1,"Internal Error: invalid 'type' (%d) in SequenceWriter::SequenceWriter()", type);
		}
	}
public:
	void write ( const FASTX& fastx, unsigned int offset=0, int len=-1 )
	{
		writer->write(fastx, offset, len);
	};
};

typedef std::auto_ptr<SequenceWriter> AutoSequenceWriter;

#endif
