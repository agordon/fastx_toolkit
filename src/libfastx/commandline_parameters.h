#ifndef __COMMANDLINE_H__
#define __COMMANDLINE_H__

#include <vector>
#include <string>
#include <getopt.h>

class CommandlineParameters
{
private:
	std::string short_options;
	std::vector<struct option> long_options;

public:
	CommandlineParameters();
	virtual ~CommandlineParameters();

	void parse_command_line ( int argc, char* argv[] ) ;

	void add_option_short_and_long(
			const std::string& long_option,
			char short_option,
			int has_arg ) ;

	void add_option_short( char short_option, int has_arg ) ;

	void add_option_long(
			const std::string& long_option,
			int long_option_value,
			int has_arg ) ;

	virtual void parameter_action(const int short_option_value, const std::string& optarg)=0;
	virtual void post_action( const std::vector<std::string> &remaining_params ) = 0;

private:
	void free_long_options();
};

#endif
