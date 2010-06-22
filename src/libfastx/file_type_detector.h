#ifndef __FILE_TYPE_DETECTOR_H__
#define __FILE_TYPE_DETECTOR_H__

#include <string>
#include <vector>

bool is_printable_string(const std::string &line);
bool is_nucleotide_string(const std::string &line);
bool is_fasta_id_string(const std::string &line);
bool is_fastq_id1_string(const std::string &line);
bool is_fastq_id2_string(const std::string &line);

bool file_type_is_readable(const std::string &filename);

bool file_type_is_fasta_single_line(const std::string &filename);
bool file_type_is_fastq(const std::string &filename);
bool file_type_is_fastx(const std::string &filename);

bool file_type_is_SAM(const std::string &filename);
bool file_type_is_BAM(const std::string &filename);

void file_type_peek_lines(const std::string &filename, size_t num_lines, std::vector< std::string >& /*output*/ lines);

#endif
