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

	virtual std::string filename() const {return _filename;}
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



/*
   Paired-end reader/writer
 */

class PE_FastqFileReader : public ISequenceReaderPE
{
	FastqFileReader end1;
	FastqFileReader end2;

	int _ASCII_quality_offset ;

public:
	PE_FastqFileReader ( const std::string& filename1, const std::string& filename2, int ASCII_quality_offset );
	PE_FastqFileReader ( input_stream_wrapper file1, input_stream_wrapper file2, int ASCII_quality_offset );

	virtual bool read_next_sequence(Sequence& /*output*/, Sequence& /*output*/) ;

	virtual ISequenceWriterPE* create_fastx_writer(const std::string& filename1, const std::string &filename2) ;
	virtual ISequenceWriterPE* create_tabular_writer(const std::string& filename) ;
};


class PE_FastqFileWriter: public ISequenceWriterPE
{
	FastqFileWriter end1;
	FastqFileWriter end2;
public:
	PE_FastqFileWriter(const std::string &filename1, const std::string &filename2);

	void write_sequence(const Sequence&, const Sequence&) ;
};
#endif
