#ifndef __FASTA_FILE_H__
#define __FASTA_FILE_H__

#include <string>

class FastaFileReader : public ISequenceReader
{
private:
	std::string _filename;
	generic_input_stream input_stream;
	size_t line_number;

public:
	FastaFileReader ( const std::string& filename ) ;
	virtual bool read_next_sequence(Sequence& /*output*/) ;
};

class FastaFileWriter : public ISequenceWriter
{
private:
	std::string _filename;
	generic_output_stream output_stream;
public:
	FastaFileWriter ( const std::string& filename ) ;
	virtual void write_sequence(const Sequence&) ;
};

#endif
