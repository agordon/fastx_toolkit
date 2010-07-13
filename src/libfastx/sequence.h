#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <vector>
#include <string>

class Sequence
{
public:
	std::string id;
	std::string nucleotides;
	std::string id2;

	std::string quality_cached_line;

	std::vector<int> quality;

	int ASCII_quality_offset;
	bool ASCII_quality_scores ; //true=ASCII, false=numeric

	Sequence();
	virtual ~Sequence() { }

	void parse_quality_scores()
	{
	}

	void clear();

	void convert_numeric_quality_score_line ( const std::string &numeric_quality_line, const std::string& _filename, size_t line_nmber);
	void convert_ascii_quality_score_line ( const std::string& numeric_quality_line );
};

/*
   Single-End Reader/Writer
*/
class ISequenceWriter;

class ISequenceReader
{
public:
	virtual ~ISequenceReader() { } ;

	virtual bool read_next_sequence(Sequence& /*output*/) = 0 ;
	virtual ISequenceWriter* create_fastx_writer(const std::string& filename) = 0 ;
	virtual ISequenceWriter* create_tabular_writer(const std::string& filename) = 0 ;

	virtual std::string filename() const = 0;
};

class ISequenceWriter
{
public:
	virtual ~ISequenceWriter() { } ;

	virtual void write_sequence(const Sequence&) = 0 ;
};



/*
   Paired-End Reader/Writer
*/
class ISequenceWriterPE;

class ISequenceReaderPE
{
public:
	virtual ~ISequenceReaderPE() { } ;

	virtual bool read_next_sequence(Sequence& /*output*/, Sequence& /*output*/) = 0 ;

	virtual ISequenceWriterPE* create_fastx_writer(const std::string& filename1, const std::string& filename2) = 0 ;
	virtual ISequenceWriterPE* create_tabular_writer(const std::string& filename) = 0 ;
};

class ISequenceWriterPE
{
public:
	virtual ~ISequenceWriterPE() { } ;

	virtual void write_sequence(const Sequence&, const Sequence&) = 0 ;
};

#endif
