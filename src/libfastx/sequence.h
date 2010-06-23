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
	std::vector<int> quality;

	int ASCII_quality_offset;
	bool ASCII_quality_scores ; //true=ASCII, false=numeric

	Sequence();
	virtual ~Sequence() { }

	void clear();
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
	virtual ISequenceWriterPE* create_writer(const std::string& filename) = 0 ;
};

class ISequenceWriterPE
{
public:
	virtual ~ISequenceWriterPE() { } ;

	virtual void write_sequence(const Sequence&, const Sequence&) = 0 ;
};

#endif
