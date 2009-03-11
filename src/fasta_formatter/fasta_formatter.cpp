#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

#include <gtextutils/stream_wrapper.h>
#include <gtextutils/text_line_reader.h>
#include <gtextutils/print_utils.h>

#include "sequence_writers.h"

using namespace std;

string input_filename;
string output_filename;
bool flag_output_empty_sequences = true;


int main()
{
	ios::sync_with_stdio(false);

	InputStreamWrapper input ( input_filename ) ;
	OutputStreamWrapper output ( output_filename );
	TextLineReader reader ( input.stream() ) ;

	SequencesWriter * pWriter = NULL ;
	SingleLineFastaWriter* pSingleLineWriter = new SingleLineFastaWriter( output.stream() ) ;
	pWriter = pSingleLineWriter;

	MultiLineFastaWriter *pMultiLineWriter = new MultiLineFastaWriter ( output.stream(), 55 ) ;
	pWriter = pMultiLineWriter ;


	if (!flag_output_empty_sequences) {
		EmptySequencesFilter *filter = new EmptySequencesFilter ( *pWriter ) ;
		pWriter = filter ;
	}
	
	int max_length = 0 ;
	string sequence_id ;
	string sequence_bases ;
	bool first_line = true ;
	while ( reader.next_line() ) {

		const string &line = reader.line_string();
		
		if ( line.length()==0 )
			continue;

		if ( line[0] == '>' ) {
			//Got new sequence identifier - print previous sequence
			if (first_line)
				first_line = false;
			else
				pWriter->write ( sequence_id, sequence_bases ) ;
			
			// Start new sequence 
			sequence_id = line ;
			sequence_bases.clear();
			sequence_bases.resize ( max_length * 2 ) ;
		} else {
			//Got sequence nucleotides
			sequence_bases += line ;
		}
	}

	//Write the last sequence
	pWriter->write ( sequence_id, sequence_bases ) ;
}

