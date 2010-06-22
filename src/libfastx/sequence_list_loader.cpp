#include <err.h>
#include <string>
#include <iterator>

#include <gtextutils/generic_input_stream.h>
#include <gtextutils/string_tokenize.h>

#include "file_type_detector.h"
#include "sequence_list_loader.h"

using namespace std;
using namespace std::tr1;

void load_sequence_ids_fasta(const std::string& filename, ISequenceIDContainer *pContainer)
{
	size_t line_number = 1 ;
	generic_input_stream input(filename);
	string id;
	while (getline(input, id)) {
		string nuc;
		if (!is_fasta_id_string(id)) {
			cerr << "Input error: Invalid FASTA ID value (" << id << ") in '" << filename
				<< "' line " << line_number << endl;
			exit(1);
		}
		const string &id_no_prefix ( id.substr(1) ) ;
		pContainer->add_sequence_id(id_no_prefix);
		++line_number;

		if (!getline(input, nuc)) {
			cerr << "Input error: failed to read nucleotides from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}

		if (!is_nucleotide_string(nuc)) {
			cerr << "Input error: Invalid nucleotides line (" << nuc << ") in '" << filename
				<< "' line " << line_number << endl;
			exit(1);
		}
		++line_number;

	}
}

void load_sequence_ids_fastq(const std::string& filename, ISequenceIDContainer *pContainer)
{
	size_t line_number = 1 ;
	generic_input_stream input(filename);
	while (input) {
		string id;
		string nuc;
		string id2;
		string quality;

		//Line 1 - ID
		if (!getline(input, id)) {
			if (input.eof())
				break;

			cerr << "Input error: failed to read ID from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}
		if (!is_fastq_id1_string(id)) {
			cerr << "Input error: Invalid FASTQ ID value (" << id << ") in '" << filename
				<< "' line " << line_number << endl;
			exit(1);
		}
		const string &id_no_prefix ( id.substr(1) ) ;
		pContainer->add_sequence_id(id_no_prefix);
		++line_number;

		//line 2 - nucleotides
		if (!getline(input, nuc)) {
			cerr << "Input error: failed to read nucleotides from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}

		if (!is_nucleotide_string(nuc)) {
			cerr << "Input error: Invalid nucleotides line (" << nuc << ") in '" << filename
				<< "' line " << line_number << endl;
			exit(1);
		}
		++line_number;

		//Line 3 - ID22
		if (!getline(input, id2)) {
			cerr << "Input error: failed to read ID-2 from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}
		if (!is_fastq_id2_string(id2)) {
			cerr << "Input error: Invalid FASTQ ID-2 value (" << id << ") in '" << filename
				<< "' line " << line_number << endl;
			exit(1);
		}
		++line_number;

		//Line 4 - quality scores (TODO: add validation)
		if (!getline(input, quality)) {
			cerr << "Input error: failed to read quality-scores from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}
		++line_number;
	}
}

void load_sequence_ids_txt(const std::string& filename, ISequenceIDContainer *pContainer)
{
	size_t line_number = 1 ;
	generic_input_stream input(filename);
	while (input) {
		string line;

		if (!getline(input, line)) {
			if (input.eof())
				break;
			cerr << "Input error: failed to read line from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}

		pContainer->add_sequence_id(line);
		++line_number;
	}
}

void load_sequence_ids_tabular(const std::string& filename, size_t column, ISequenceIDContainer *pContainer)
{
	size_t line_number = 1 ;
	generic_input_stream input(filename);
	while (input) {
		string line;

		if (!getline(input, line)) {
			if (input.eof())
				break;

			cerr << "Input error: failed to read line from '" << filename
				<< "' line " << line_number << ":" << string_error(errno) << endl;
			exit(1);
		}

		if (line.size()==0 || line.at(0)=='#') {
			++line_number;
			continue;
		}

		vector<string> fields;
		String_Tokenize ( line, std::back_inserter(fields), "\t" );

		if (fields.size()<=column) {
			cerr << "Input error: Not enough fields in '" << filename
				<< "' line " << line_number << ": expecting " << column+1 << " got only " << fields.size() << endl;
			exit(1);
		}

		const string value(fields[column]);
		pContainer->add_sequence_id(value);

		++line_number;
	}
}

void load_sequence_ids(LISTTYPE input_type, const std::string& filename, size_t optional_column, ISequenceIDContainer *pContainer)
{
	if (input_type==TYPE_AUTO_DETECT) {
		if (file_type_is_fastx(filename))
			input_type=TYPE_FASTX;
/*		else if (file_type_is_SAM(filename))
			input_type=TYPE_SAM;
		else if (file_type_is_BAM(filename))
			input_type=TYPE_BAM; */
		else
			input_type=TYPE_SIMPLE_TEXT;
	}

	switch(input_type)
	{
	case TYPE_FASTX:
		if (file_type_is_fasta_single_line(filename))
			load_sequence_ids_fasta(filename,pContainer);
		else
		if (file_type_is_fastq(filename))
			load_sequence_ids_fastq(filename,pContainer);
		else
			errx(1,"Internal error: input_type=TYPE_FASTX but file is not detected as FASTQ or FASTA");
		break;

	case TYPE_SIMPLE_TEXT:
		load_sequence_ids_txt(filename,pContainer);
		break;

	case TYPE_TABULAR:
		load_sequence_ids_tabular(filename, optional_column, pContainer);
		break;

	case TYPE_SAM:
		errx(1,"Internal error: SAM input is not yet implemented");
		break;

	case TYPE_BAM:
		errx(1,"Internal error: BAM input is not yet implemented");
		break;

	case TYPE_AUTO_DETECT:
		errx(1,"Internal error: input_type=AUTO_DETECT at " __FILE__ ": %d", __LINE__);
		break;
	default:
		errx(1,"Internal error: invalid input_type='%d'", (int)input_type);
		break;
	}
}

