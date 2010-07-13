#ifndef __FASTA_FILE_H__
#define __FASTA_FILE_H__

#include <string>

/*
   Single-End Reader/Writer
 */

class FastaFileReader : public ISequenceReader
{
private:
	std::string _filename;
	generic_input_stream input_stream;
	size_t line_number;

public:
	FastaFileReader ( const std::string& filename ) ;
	FastaFileReader ( input_stream_wrapper w ) ;
	virtual bool read_next_sequence(Sequence& /*output*/) ;

	virtual ISequenceWriter* create_fastx_writer(const std::string& filename);
	virtual ISequenceWriter* create_tabular_writer(const std::string& filename);
	virtual std::string filename() const {return _filename;}
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


/*
   Paired-end reader/writer
 */

class PE_FastaFileReader : public ISequenceReaderPE
{
	FastaFileReader end1;
	FastaFileReader end2;

public:
	PE_FastaFileReader ( const std::string& filename1, const std::string& filename2 );
	PE_FastaFileReader ( input_stream_wrapper file1, input_stream_wrapper file2 );

	virtual bool read_next_sequence(Sequence& /*output*/, Sequence& /*output*/) ;

	virtual ISequenceWriterPE* create_fastx_writer(const std::string& filename1, const std::string &filename2) ;
	virtual ISequenceWriterPE* create_tabular_writer(const std::string& filename) ;
};


class PE_FastaFileWriter: public ISequenceWriterPE
{
	FastaFileWriter end1;
	FastaFileWriter end2;
public:
	PE_FastaFileWriter(const std::string &filename1, const std::string &filename2);

	void write_sequence(const Sequence&, const Sequence&) ;
};

#endif
