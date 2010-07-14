#include <climits>
#include <string>
#include <iostream>

#include "libfastx/fastxse_commandline_parameters.h"

using namespace std;

class FastxCollapserCommandLine : public FastxSE_commandline_parameters
{
private:
public:
	void print_help()
	{
		cout << "No help for you, come back - one year!" << endl;
	}
};

int main(int argc, char* argv[])
{
	FastxCollapserCommandLine p;
	p.parse_command_line(argc,argv);

	Sequence seq;
	ISequenceReader *reader = p.reader().get();
	ISequenceWriter *writer = p.writer().get();

	while (reader->read_next_sequence(seq)) {
		writer->write_sequence(seq);
	}

	return 0;
}
