#ifndef __FASTXSE_COMMANDLINE_PARAMETERS__
#define __FASTXSE_COMMANDLINE_PARAMETERS__

#include <climits>
#include <tr1/memory>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>

#include "commandline_parameters.h"
#include "sequence.h"
#include "tab_file.h"
#include "fastx_file.h"

enum {
	OPT_TAB_IN = CHAR_MAX+6789,
	OPT_TAB_OUT
};

class ISequenceReader;
class ISequenceWriter;

typedef std::tr1::shared_ptr<ISequenceReader> SharedSequenceReader;
typedef std::tr1::shared_ptr<ISequenceWriter> SharedSequenceWriter;

class FastxSE_commandline_parameters : public CommandlineParameters
{
private:
	std::string _input_filename;
	std::string _output_filename;
	int	    _verbose;
	bool	    _tabular_input;
	bool	    _tabular_output;
	int	    _ASCII_quality_offset;

	SharedSequenceReader _Reader;
	SharedSequenceWriter _Writer;

public:
	FastxSE_commandline_parameters();

	void parameter_action(const int short_option_value, const std::string& optarg);
	void post_action( const std::vector<std::string> &remaining_params );

	virtual void print_help()=0;

	virtual int verbose() const { return _verbose; }

	SharedSequenceReader reader()
	{
		if (_Reader.get()==NULL) {
			if (_tabular_input) {
				_Reader = SharedSequenceReader(new TabularFileReader(_input_filename, _ASCII_quality_offset));
			} else {
				_Reader = SharedSequenceReader(create_fastx_reader(_input_filename, _ASCII_quality_offset));
			}
		}
		return _Reader;
	}

	SharedSequenceWriter writer()
	{
		if (_Writer.get()==NULL) {
			if (_tabular_output) {
				_Writer = SharedSequenceWriter ( reader()->create_tabular_writer(_output_filename) );
			} else {
				_Writer = SharedSequenceWriter ( reader()->create_fastx_writer(_output_filename) ) ;
			}
		}
		return _Writer;
	}

	std::ostream& verbose_stream() const
	{
		if (_output_filename.empty())
			//output goes to STDOUT, verbose report to STDERR
			return std::cerr;
		else
			//output goes to a file, vebose report goes to STDOUT
			return std::cout;
	}
};


#endif
