#include <climits>
#include <string>
#include <iostream>

#include "libfastx/commandline_parameters.h"

using namespace std;

class params : public CommandlineParameters
{
public:
	virtual void parameter_action(const int short_option_value, const std::string& optarg)
	{
		if (isprint(short_option_value)) {
			cout << "option '" << (char)short_option_value << "'" ;
		} else
		if (short_option_value > CHAR_MAX) {
			cout << "option CHAR_MAX+" << (int)(short_option_value-CHAR_MAX);
		} else {
			cout << "option " << short_option_value << "(warning: strange value?)" ;
		}

		if (!optarg.empty()) {
			cout << "\targ = '" << optarg << "'" ;
		}
		cout << endl;
	}

	virtual void post_action( const std::vector<std::string> &remaining_params )
	{
		std::vector<std::string>::const_iterator it = remaining_params.begin();
		while (  it != remaining_params.end() ) {
			cout << "Extra parameters: " << (*it) << endl;
			++it;
		}
	}
};

int main(int argc, char* argv[])
{
	params p;

	p.add_option_short_and_long("help", 'h', no_argument);
	p.add_option_short_and_long("verbose", 'v', no_argument);
	p.add_option_short_and_long("input", 'i', required_argument ) ;
	p.add_option_short('l', required_argument);
	p.add_option_short('d', no_argument);
	p.add_option_long("qual", CHAR_MAX+5, no_argument);
	p.parse_command_line(argc, argv);

	return 0;
}
