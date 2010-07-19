#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <vector>
#include <string>
#include <sstream>

class Sequence
{
public:
	std::string id;
	std::string nucleotides;
	std::string id2;

	std::string quality_cached_line;

	//std::vector<int> quality;
	int ASCII_quality_offset;
	bool ASCII_quality_scores ; //true=ASCII, false=numeric
	mutable int cached_multiplicity_count;

	Sequence();
	virtual ~Sequence() { }

	void parse_quality_scores()
	{
	}

	size_t get_multiplicity_count() const
	{
		if (cached_multiplicity_count!=0)
			return cached_multiplicity_count;

		size_t i = id.find_last_of('-');
		//can't detect multiplicity ? assume it's one.
		if (i==std::string::npos) {
			cached_multiplicity_count =1 ;
			return 1;
		}
		std::stringstream ss(id.substr(i+1));
		size_t count;
		ss >> count;
		if (ss.good()) {
			cached_multiplicity_count = count;
			return count;
		}

		cached_multiplicity_count = 1 ;
		return 1;
	}

	/*
	   returns:
	      1 - if this is END1 of a paired-end (i.e. "ID" string ends in "/1")
	      2 - if this is END2 of a paired-end (i.e. "ID" string ends in "/2")
	      0 - could not detect END, probably not part of a paired-end FASTQ

	 */
	size_t get_paired_end()
	{
		if (id.length()<3)
			return 0;
		if (id[id.length()-2] != '/')
			return 0;
		const char c = id[id.length()-1];
		if (c=='1')
			return 1;
		if (c=='2')
			return 2;
		return 0;
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
