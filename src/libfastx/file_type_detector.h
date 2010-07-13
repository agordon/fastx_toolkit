#ifndef __FILE_TYPE_DETECTOR_H__
#define __FILE_TYPE_DETECTOR_H__

#include <string>
#include <vector>

#define SKIP_PRINTABLE_STRING_CHECK
#define INLINE_STRING_FUNCTIONS

#ifdef INLINE_STRING_FUNCTIONS

inline bool is_printable_string(const std::string &
#ifndef SKIP_PRINTABLE_STRING_CHECK
		line
#endif
		)
{
#ifndef SKIP_PRINTABLE_STRING_CHECK
	for (size_t i=0;i<line.size();++i) {
		const char c = line[i];
		if ( ! (isprint(c) || isblank(c)) )
			return false;
	}
#endif
	return true;
}

inline bool is_nucleotide_string(const std::string &line)
{
	if (line.empty())
		return false;

	size_t count = line.length();
	const char *pC = line.data();
	while ( count ) {
		switch(*pC)
		{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'N':
			break;

		default:
			return false;
		}
		++pC;
		--count;
	}
	return true;
}

inline bool is_fasta_id_string(const std::string &line)
{
	if (line.size()<2)
		return false;
	if (line[0] != '>')
		return false;
	return is_printable_string(line);
}

inline bool is_fastq_id1_string(const std::string &line)
{
	if (line.empty())
		return false;
	if (line[0] != '@')
		return false;
	return is_printable_string(line);
}

inline bool is_fastq_id2_string(const std::string &line)
{
	if (line.empty())
		return false;
	if (line[0] != '+')
		return false;
	if (line.size()==1) //on some machines, the second ID in a FASTQ file can be empty
		return true;
	return is_printable_string(line);
}
#else
bool is_printable_string(const std::string &line);
bool is_nucleotide_string(const std::string &line);
bool is_fasta_id_string(const std::string &line);
bool is_fastq_id1_string(const std::string &line);
bool is_fastq_id2_string(const std::string &line);
#endif


bool file_type_is_readable(const std::string &filename);
bool file_type_is_fasta_single_line(const std::string &filename);
bool file_type_is_fastq(const std::string &filename);
bool file_type_is_fastx(const std::string &filename);

bool file_type_is_SAM(const std::string &filename);
bool file_type_is_BAM(const std::string &filename);

void file_type_peek_lines(const std::string &filename, size_t num_lines, std::vector< std::string >& /*output*/ lines);

#endif
