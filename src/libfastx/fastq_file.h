#ifndef __FASTQ_FILE_H__
#define __FASTQ_FILE_H__

#include <string>

class FastqFileReader : public ISequenceReader
{
private:
	std::string _filename;
	generic_input_stream input_stream;
	size_t line_number;

	int _ASCII_quality_offset ;

public:
	FastqFileReader ( const std::string& filename , int ASCII_quality_offset ) ;
	FastqFileReader ( input_stream_wrapper w, int ASCII_quality_offset ) ;
	virtual bool read_next_sequence(Sequence& /*output*/) ;

	virtual ISequenceWriter* create_fastx_writer(const std::string& filename);
	virtual ISequenceWriter* create_tabular_writer(const std::string& filename);

public:
	void convert_numeric_quality_score_line ( const std::string& numeric_quality_line, std::vector<int>& /*output*/);
	static void convert_ascii_quality_score_line ( const std::string& numeric_quality_line, std::vector<int>& /*output*/, int ASCII_OFFSET);
};

class FastqFileWriter : public ISequenceWriter
{
private:
	std::string _filename;
	generic_output_stream output_stream ;
public:
	FastqFileWriter ( const std::string& filename) ;
	virtual void write_sequence(const Sequence&) ;

};

#endif
