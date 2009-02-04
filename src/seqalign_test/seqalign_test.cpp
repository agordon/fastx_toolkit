#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include "sequence_alignment.h"


int main( /*int argc, char* argv[] */)
{
	HalfLocalSequenceAlignment lsa ;

	const SequenceAlignmentResults& results = lsa.align("AAAGGTTTCCC","AGGCTT" );
	lsa.print_matrix();
	results.print();


	return 0;
}
