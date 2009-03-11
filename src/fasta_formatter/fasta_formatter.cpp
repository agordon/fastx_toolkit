#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

#include <gtextutils/stream_wrapper.h>
#include <gtextutils/text_line_reader.h>
#include <gtextutils/print_utils.h>

using namespace std;

string input_filename;
string output_filename;
bool flag_output_empty_sequences = true;

class SequencesWriter
{
public:
	virtual ~SequencesWriter() {}
	virtual void write ( const std::string& sequence_id, const std::string& sequence_bases ) = 0 ;
};

class EmptySequencesFilter : public SequencesWriter
{
private:
	SequencesWriter& upstream ;

public:
	EmptySequencesFilter ( SequencesWriter & _upstream ) : upstream(_upstream) {}

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases)
	{
		if ( !sequence_bases.empty() )
			upstream.write ( sequence_id, sequence_bases ) ;
	}
};

class SingleLineFastaWriter : public SequencesWriter
{
private:
	ostream& ostrm ;
public:
	SingleLineFastaWriter ( ostream& output_stream ) :
		ostrm ( output_stream )
	{
	}

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id << endl;
		if ( !sequence_bases.empty() )
			ostrm << sequence_bases << endl ;
	}
};

class MultiLineFastaWriter : public SequencesWriter
{
private:
	ostream& ostrm ;
	size_t   max_width ;

public:
	MultiLineFastaWriter ( ostream& output_stream, size_t _max_width ) :
		ostrm ( output_stream ), max_width ( _max_width ) 
	{
	}

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id << endl;
		if ( !sequence_bases.empty() ) {
			size_t start = 0 ;
			while ( (sequence_bases.length() - start) > max_width ) {
				ostrm << sequence_bases.substr ( start, max_width ) << endl;
				start += max_width ;
			}
			if ( sequence_bases.length() - start > 0 )
				ostrm << sequence_bases.substr ( start ) << endl ;
		}
	}
};

class TabulatedFastaWriter : public SequencesWriter
{
private:
	ostream& ostrm ;
public:
	TabulatedFastaWriter ( ostream& output_stream ) : ostrm ( output_stream ) { }

	virtual void write ( const std::string & sequence_id, const std::string& sequence_bases )
	{
		ostrm << sequence_id ;
		if ( !sequence_bases.empty() ) {
			ostrm << "\t" ;
			ostrm << sequence_bases ;
		}
		ostrm << endl;
	}
};

int main()
{
	ios::sync_with_stdio(false);

	InputStreamWrapper input ( input_filename ) ;
	OutputStreamWrapper output ( output_filename );
	TextLineReader reader ( input.stream() ) ;

	SequencesWriter * pWriter = NULL ;
	SingleLineFastaWriter* pSingleLineWriter = new SingleLineFastaWriter( output.stream() ) ;
	pWriter = pSingleLineWriter;

	MultiLineFastaWriter *pMultiLineWriter = new MultiLineFastaWriter ( output.stream(), 55 ) ;
	pWriter = pMultiLineWriter ;


	if (!flag_output_empty_sequences) {
		EmptySequencesFilter *filter = new EmptySequencesFilter ( *pWriter ) ;
		pWriter = filter ;
	}
	
	int max_length = 0 ;
	string sequence_id ;
	string sequence_bases ;
	bool first_line = true ;
	while ( reader.next_line() ) {

		const string &line = reader.line_string();
		
		if ( line.length()==0 )
			continue;

		if ( line[0] == '>' ) {
			//Got new sequence identifier - print previous sequence
			if (first_line)
				first_line = false;
			else
				pWriter->write ( sequence_id, sequence_bases ) ;
			
			// Start new sequence 
			sequence_id = line ;
			sequence_bases.clear();
			sequence_bases.resize ( max_length * 2 ) ;
		} else {
			//Got sequence nucleotides
			sequence_bases += line ;
		}
	}

	//Write the last sequence
	pWriter->write ( sequence_id, sequence_bases ) ;
}

