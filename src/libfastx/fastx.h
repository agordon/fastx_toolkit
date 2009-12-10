/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __FASTX_HEADER__
#define __FASTX_HEADER__

#ifdef __cplusplus
extern "C" {
#endif

/* for PATH_MAX */
#include <limits.h>

#define MIN_QUALITY_VALUE (-50)
#define MAX_QUALITY_VALUE 50
#define QUALITY_VALUES_RANGE (MAX_QUALITY_VALUE-MIN_QUALITY_VALUE)


#ifndef MAX_SEQ_LINE_LENGTH 
#define MAX_SEQ_LINE_LENGTH (25000)
#endif

typedef enum {
	FASTA_ONLY=0,
	FASTA_OR_FASTQ=1,
	FASTQ_ONLY=2	
} ALLOWED_INPUT_FILE_TYPES;

typedef enum {
	DISALLOW_N=0,
	ALLOW_N=1,
	ALLOW_U=2	
} ALLOWED_INPUT_BASES;

typedef enum {
	REQUIRE_UPPERCASE=0,
	ALLOW_LOWERCASE=1
} ALLOWED_INPUT_CASE;

typedef enum {
	OUTPUT_FASTA=0,
	OUTPUT_FASTQ_ASCII_QUAL=1,
	OUTPUT_FASTQ_NUMERIC_QUAL=2,
	OUTPUT_SAME_AS_INPUT=3
} OUTPUT_FILE_TYPE;

#pragma pack(1) 
typedef struct 
{
	/* Record data - common for FASTA/FASTQ */
	char    input_sequence_id_prefix[1];   //DON'T touch this - this hack will read the entire name into the variable 'name',
				  //leaving the prefix ('>' or '@') in 'input_sequence_id_name'.
	char    name[MAX_SEQ_LINE_LENGTH+1];
	char    nucleotides[MAX_SEQ_LINE_LENGTH+1];
	/* Record data - only for FASTQ */
	char    input_name2_prefix[1];         //same hack as 'input_sequence_id_prefix'
	char	name2[MAX_SEQ_LINE_LENGTH+1];
	int	quality[MAX_SEQ_LINE_LENGTH+1];  //note: this is NOT ascii values, but numerical values
					       //      numeric quality scores and ASCII quality scores
					       //      are automatically converted to numbers (-15 to 93)

	/* Configuration */
	int	allow_input_filetype;	// 0 = Allow only FASTA
	int	allow_N;		// 1 = N is valid nucleotide, 0 = only A/G/C/T are valid
	int	allow_U;		// 1 = U is a valid nucleotide
	int	allow_lowercase;	
	int	read_fastq;		// 1 = Input is FASTQ (only if allow_input_fastq==1)
	int	read_fastq_ascii;	// 1 = Input is FASTQ with ASCII quality scores (0 = with numeric quality scores)
	int	write_fastq;		// 0 = Write only FASTA (regardless of input type)
	int	write_fastq_ascii;	// 1 = Write ASCII quality scores, 0 = write numeric quality scores
	int	compress_output;		// 1 = pass output through GZIP

	int     copy_input_fastq_format_to_output ; // 1 = copy 'read_fastq_ascii' to 'write_fastq_ascii'
						    // so that the output format is the same as the input
	int	fastq_ascii_quality_offset ;  // When reading FASTQ with ASCII quality values, subtract this value to get 'real' numeric quality score;


	/* Internal data */
	int	allowed_nucleotides[256];	//quick lookup table for valid input	
	char	output_sequence_id_prefix;	// '>' or '@', depending on the requested output type

	char	input_file_name[PATH_MAX];	//in linux, PATH_MAX is defined in <linux/limits.h>
	unsigned long long input_line_number;
	char	output_file_name[PATH_MAX];	//in linux, PATH_MAX is defined in <linux/limits.h>

	size_t	num_input_sequences;
	size_t  num_output_sequences;
	size_t  num_input_reads;
	size_t  num_output_reads;

	FILE*	input;
	FILE*	output;
} FASTX ;


void fastx_init_reader(FASTX *pFASTX, const char* filename, 
		ALLOWED_INPUT_FILE_TYPES allowed_input_filetype,
		ALLOWED_INPUT_BASES allow_bases,
		ALLOWED_INPUT_CASE allow_lowercase,
		int fastq_ascii_quality_offset);

// If the sequence identifier is collapsed (= "N-N") returns the reads_count,
// otherwise, returns 1
int get_reads_count(const FASTX *pFASTX);

void fastx_init_writer(FASTX *pFASTX,
		const char* filename,
		OUTPUT_FILE_TYPE output_type,
		int compress_output);
	
int fastx_read_next_record(FASTX *pFASTX);

void fastx_write_record(FASTX *pFASTX);

size_t num_input_sequences(const FASTX *pFASTX);
size_t num_input_reads(const FASTX *pFASTX);
size_t num_output_sequences(const FASTX *pFASTX);
size_t num_output_reads(const FASTX *pFASTX);

#ifdef __cplusplus
}
#endif

#endif

