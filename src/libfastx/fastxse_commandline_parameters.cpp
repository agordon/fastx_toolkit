#include <err.h>
#include <cstdlib>
#include "config.h"
#include "commandline_parameters.h"
#include "fastxse_commandline_parameters.h"

FastxSE_commandline_parameters::FastxSE_commandline_parameters() :
	_verbose(0), _tabular_input(false), _tabular_output(false)
{
	add_option_short_and_long("help", 'h', no_argument);
	add_option_short_and_long("verbose", 'v', no_argument);
	add_option_short_and_long("input", 'i', required_argument ) ;
	add_option_short_and_long("output", 'o', required_argument ) ;
	add_option_long("tabin", OPT_TAB_IN, no_argument ) ;
	add_option_long("tabout", OPT_TAB_OUT, no_argument ) ;
	add_option_long("pipein", OPT_TAB_IN, no_argument ) ;
	add_option_long("pipeout", OPT_TAB_OUT, no_argument ) ;
	add_option_short('Q',required_argument);
}

void FastxSE_commandline_parameters::parameter_action(const int short_option_value, const std::string& optarg)
{
	switch(short_option_value)
	{
	case 'h':
		print_help();
		exit(0);
		break;

	case 'v':
		++_verbose;
		break;

	case 'i':
		_input_filename = optarg ;
		break;

	case 'o':
		_output_filename = optarg ;
		break;

	case 'Q':
		_ASCII_quality_offset = atoi(optarg.c_str());
		if (_ASCII_quality_offset==0)
			errx(1,"Error: invalid ASCII quality offset '%s'", optarg.c_str());
		break;

	case OPT_TAB_IN:
		_tabular_input = true ;
		break;

	case OPT_TAB_OUT:
		_tabular_output = true ;
		break;

	default:
		errx(1,"Error: unknown command line option '%c'", short_option_value );
	}
}

void FastxSE_commandline_parameters::post_action( const std::vector<std::string> &remaining_params )
{
	if (_input_filename.empty() && !remaining_params.empty()) {
		_input_filename = remaining_params[0];
	}
}

