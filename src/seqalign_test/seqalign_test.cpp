#include <string>
#include <vector>
#include "sequence_alignment.h"


int main(int argc, char* argv[])
{
	LocalSequenceAlignment lsa ;

	const SequenceAlignmentResults& results = lsa.align("AAAGGTTTCCC","AGGCTT" );
	lsa.print_matrix();
	results.print();


	return 0;
}
