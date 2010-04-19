#ifndef __FASTX_FILE_WRITER__
#define __FASTX_FILE_WRITER__

#include <ostream>
#include <string>

#include <gtextutils/stream_wrapper.h>

#include "fastx.h"

class FastxWriter
{
public:
	virtual void write ( const FASTX& fastx, unsigned int offset, int len )=0;
};

class FastqWriter : public FastxWriter
{
private:
	std::ostream &out;
	std::string name_suffix;

public:
	FastqWriter ( std::ostream& output_stream, const std::string &_name_suffix );
	void write ( const FASTX& fastx, unsigned int offset, int len );
};

class FastaWriter : public FastxWriter
{
private:
	std::ostream &out;
	std::string name_suffix;

public:
	FastaWriter ( std::ostream& output_stream, const std::string &_name_suffix );
	void write ( const FASTX& fastx, unsigned int offset, int len );
};

#endif
