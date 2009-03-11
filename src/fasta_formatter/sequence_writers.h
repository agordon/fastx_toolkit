#ifndef __SEQUENCE_WRITERS__
#define __SEQUENCE_WRITERS__

#include <string>
#include <ostream>

class SequencesWriter
{
public:
	virtual ~SequencesWriter() {}
	virtual void write ( const std::string& sequence_id, const std::string& sequence_bases ) = 0 ;
};

class EmptySequencesFilter : public SequencesWriter
{
private:
	SequencesWriter* upstream ;

public:
	EmptySequencesFilter ( SequencesWriter * _upstream ) : upstream(_upstream) {}

	~EmptySequencesFilter()
	{
		delete upstream;
	}

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases)
	{
		if ( !sequence_bases.empty() )
			upstream->write ( sequence_id, sequence_bases ) ;
	}
};

class SingleLineFastaWriter : public SequencesWriter
{
private:
	std::ostream& ostrm ;
public:
	SingleLineFastaWriter ( std::ostream& output_stream ) : ostrm ( output_stream ) { }

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id << std::endl;
		if ( !sequence_bases.empty() )
			ostrm << sequence_bases << std::endl ;
	}
};

class MultiLineFastaWriter : public SequencesWriter
{
private:
	std::ostream& ostrm ;
	size_t   max_width ;

public:
	MultiLineFastaWriter ( std::ostream& output_stream, size_t _max_width ) :
		ostrm ( output_stream ), max_width ( _max_width ) 
	{
	}

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id << std::endl;
		if ( !sequence_bases.empty() ) {
			size_t start = 0 ;
			while ( (sequence_bases.length() - start) >= max_width ) {
				ostrm << sequence_bases.substr ( start, max_width ) << std::endl;
				start += max_width ;
			}
			if ( sequence_bases.length() - start > 0 )
				ostrm << sequence_bases.substr ( start ) << std::endl ;
		}
	}
};

class TabulatedFastaWriter : public SequencesWriter
{
private:
	std::ostream& ostrm ;
public:
	TabulatedFastaWriter ( std::ostream& output_stream ) : ostrm ( output_stream ) { }

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id.substr(1) ;
		if ( !sequence_bases.empty() ) {
			ostrm << "\t" ;
			ostrm << sequence_bases ;
		}
		ostrm << std::endl;
	}
};

#endif

