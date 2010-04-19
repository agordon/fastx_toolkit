#include <err.h>
#include "fastx_writer.h"

using namespace std;

FastqWriter::FastqWriter ( std::ostream& output_stream , const std::string &_name_suffix ) :
	out ( output_stream ),
	name_suffix ( _name_suffix )
{
}

inline
std::string trim_string(const std::string& input, unsigned int offset, int len)
{
	//Todo: make this into an exception
	if (len==0)
		errx(1, "Internal error: len==0 in trim_strin(input='%s', offset=%u)",
				input.c_str(), offset);


	if (offset==0 && len==-1)
		return input;

	if (offset>0 && len==-1)
		return input.substr(offset);

	return input.substr(offset,len);
}

void FastqWriter::write ( const FASTX& fastx, unsigned int offset, int len )
{
	out << "@" << fastx.name ;
	if ( ! name_suffix.empty() )
		out << name_suffix;
	out << endl;

	out <<  trim_string(fastx.nucleotides, offset,len) << endl;

	out << "+" << fastx.name2 ;
	if ( ! name_suffix.empty() )
		out << name_suffix;
	out << endl;

	//TODO: convert qualities to ASCII,
	//      with correct OFFSET
	out << "TODO: Qualities output" << endl;
}

FastaWriter::FastaWriter ( std::ostream& output_stream , const std::string &_name_suffix ) :
	out ( output_stream ),
	name_suffix ( _name_suffix )
{
}

void FastaWriter::write ( const FASTX& fastx, unsigned int offset, int len )
{
	out << ">" << fastx.name ;
	if ( ! name_suffix.empty() )
		out << name_suffix;
	out << endl;

	out <<  trim_string(fastx.nucleotides, offset,len) << endl;
}
