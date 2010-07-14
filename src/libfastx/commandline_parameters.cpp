#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <err.h>
#include "config.h"

#include "sequence.h"
#include "commandline_parameters.h"

#include <iostream>
using namespace std;

CommandlineParameters::CommandlineParameters()
{
}

CommandlineParameters::~CommandlineParameters()
{
	free_long_options();
}

void CommandlineParameters::parse_command_line ( int argc, char* argv[] )
{
	int c;

	struct option op = { NULL,0,0,0 };
	long_options.push_back(op);

	while ( (c=getopt_long(argc,argv,short_options.c_str(), &long_options[0], NULL )) != -1 ) {
		if (c=='?') {
			exit(1);
		}
		std::string optarg_str("");
		if (optarg)
			optarg_str = optarg ;
		parameter_action(c, optarg_str) ;
	}

	std::vector<std::string> remaining_params;

	for ( int i = optind ; i < argc ; ++i )
		remaining_params.push_back(argv[i]);

	post_action ( remaining_params ) ;
}

void CommandlineParameters::add_option_short_and_long(
		const std::string& long_option,
		char short_option,
		int has_arg  )
{
	if (short_options.find(short_option) != std::string::npos)
		errx(1,"Internal error: add_option_short_and_long: short-option '%c' is already used!", short_option);
	if (!isalnum(short_option))
		errx(1,"Internal error: add_option_short_and_long: short-option 0x%02d is not an alpha-numeric ascii character!", short_option);

	struct option op;

	op.name = strdup(long_option.c_str());
	op.has_arg = has_arg;
	op.flag = NULL ;
	op.val = short_option ;
	long_options.push_back(op);

	short_options += short_option;
	if (has_arg)
		short_options += ":";
}

void CommandlineParameters::add_option_short( char short_option, int has_arg )
{
	if (short_options.find(short_option) != std::string::npos)
		errx(1,"Internal error: add_optio_shortn: short-option '%c' is already used!", short_option);
	if (!isalnum(short_option))
		errx(1,"Internal error: add_option_short: short-option 0x%02d is not an alpha-numeric ascii character!", short_option);

	short_options += short_option;
	if (has_arg)
		short_options += ":";
}

void CommandlineParameters::add_option_long(
		const std::string& long_option,
		int long_option_value,
		int has_arg)
{
	struct option op;

	op.name = strdup(long_option.c_str());
	op.has_arg = has_arg;
	op.flag = NULL ;
	op.val = long_option_value ;
	long_options.push_back(op);
}

void CommandlineParameters::free_long_options()
{
	std::vector<struct option>::iterator it;
	for (it = long_options.begin(); it != long_options.end(); ++it) {
		free((void*)it->name);
		it->name=NULL ;
	}
}
