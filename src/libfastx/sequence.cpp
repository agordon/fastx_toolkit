#include "sequence.h"

Sequence::Sequence()
{
	clear();
}

void Sequence::clear()
{
	id.clear();
	nucleotides.clear();
	id2.clear();
	quality.clear();

	ASCII_quality_offset = 0;
	ASCII_quality_scores = true;
}
