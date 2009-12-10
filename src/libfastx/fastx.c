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
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#include "chomp.h"
#include "fastx.h"

/*
	valid_sequence_string - 
		check validity of a given sequence string.
	
	input - 
		sequence - NULL terminated string to be validated.

	Output - 
		1 (true) - The given sequence is valid - contained only A/C/G/N/T characters.
		0 (false) - The given string contained invalid characeters.

	Remark -
		sequences with unknown (N) bases are considered VALID.
*/
static int validate_nucleotides_string(int allowed_nucleotides[256], const char* seq)
{
	int match = 1 ;
	while (*seq != '\0' && match) {
		match &=  allowed_nucleotides[ (int) *seq ];
		seq++;
	}
	return match;
}

static void create_lookup_table(FASTX *pFASTX)
{
	int i;

	for (i=0; i<256; i++)
		pFASTX->allowed_nucleotides[i] = 0 ;

	pFASTX->allowed_nucleotides['A'] = 1;
	pFASTX->allowed_nucleotides['C'] = 1;
	pFASTX->allowed_nucleotides['G'] = 1;
	pFASTX->allowed_nucleotides['T'] = 1;

	if (pFASTX->allow_N)
		pFASTX->allowed_nucleotides['N'] = 1;
	if (pFASTX->allow_U)
		pFASTX->allowed_nucleotides['U'] = 1;

	if (pFASTX->allow_lowercase) {
		pFASTX->allowed_nucleotides['a'] = 1;
		pFASTX->allowed_nucleotides['c'] = 1;
		pFASTX->allowed_nucleotides['g'] = 1;
		pFASTX->allowed_nucleotides['t'] = 1;

		if (pFASTX->allow_N)
			pFASTX->allowed_nucleotides['n'] = 1;
		if (pFASTX->allow_U)
			pFASTX->allowed_nucleotides['u'] = 1;
	}
}

static void detect_input_format(FASTX *pFASTX)
{
	//Get the first character in the file,
	//and put it right back
	int c = fgetc(pFASTX->input);
	ungetc(c, pFASTX->input);
	
	switch(c) {
	case '>':	/* FASTA file */
		if ( pFASTX->allow_input_filetype==FASTQ_ONLY )
			errx(1,"input file (%s) is FASTA, but only FASTQ input is allowed.", 
				pFASTX->input_file_name);
		pFASTX->read_fastq = 0 ;
		break;

	case '@':	/* FASTQ file */
		if ( pFASTX->allow_input_filetype==FASTA_ONLY )
			errx(1,"input file (%s) is FASTQ, but only FASTA input is allowed.", 
				pFASTX->input_file_name);
		pFASTX->read_fastq = 1;	
		break;
	
	case -1:   /* EOF as first character - no input */
		errx(1, "Premature End-Of-File (filename ='%s')", pFASTX->input_file_name);
		break; 

	default:
		errx(1, "input file (%s) has unknown file format (not FASTA or FASTQ), first character = %c (%d)", 
			pFASTX->input_file_name, c,c);
	}
}

static void convert_ascii_quality_score_line(const char* ascii_quality_scores, FASTX *pFASTX)
{
	size_t i;

	if (strlen(ascii_quality_scores) != strlen(pFASTX->nucleotides))
		errx(1,"number of quality values (%zu) doesn't match number of nucleotides (%zu) on line %lld",
				strlen(ascii_quality_scores), strlen(pFASTX->nucleotides),
				pFASTX->input_line_number);

	for (i=0; i<strlen(ascii_quality_scores); i++) {
		pFASTX->quality[i] = (int) (ascii_quality_scores[i] - pFASTX->fastq_ascii_quality_offset ) ;
		if (pFASTX->quality[i] < -15 || pFASTX->quality[i] > 93) 
			errx(1, "Invalid quality score value (char '%c' ord %d quality value %d) on line %lld",
				ascii_quality_scores[i], ascii_quality_scores[i],
				pFASTX->quality[i], pFASTX->input_line_number );
	}

}

static void convert_numeric_quality_score_line ( const char* numeric_quality_line, FASTX *pFASTX )
{
	size_t index;
	const char *quality_tok;
	char *endptr;
	int quality_value;

	index=0;
	quality_tok = numeric_quality_line;
	do {
		//read the quality score as an integer value
		quality_value = strtol(quality_tok, &endptr, 10);
		if (endptr == quality_tok) 
			errx(1,"Error: invalid quality score data on line %lld (quality_tok = \"%s\"", 
				pFASTX->input_line_number ,quality_tok);

		if (quality_value > 93 || quality_value < -15)
			errx(1, "invalid quality score value (%d) in line %lld.", 
				quality_value, pFASTX->input_line_number);
		
		//convert it ASCII (as per solexa's encoding)
		pFASTX->quality[index] = quality_value; 
		index++;
		quality_tok = endptr;
	} while (quality_tok != NULL && *quality_tok!='\0') ;

	if (index != strlen(pFASTX->nucleotides)) {
		errx(1,"number of quality values (%zu) doesn't match number of nucleotides (%zu) on line %lld",
				index, strlen(pFASTX->nucleotides), pFASTX->input_line_number );
	}
}

void fastx_init_reader(FASTX *pFASTX, const char* filename, 
		ALLOWED_INPUT_FILE_TYPES allowed_input_filetype,
		ALLOWED_INPUT_BASES allow_bases,
		ALLOWED_INPUT_CASE allow_lowercase,
		int fastq_ascii_quality_offset)
{
	if (pFASTX==NULL)
		errx(1,"Internal error: pFASTX==NULL (%s:%d)", __FILE__,__LINE__);

	memset(pFASTX, 0, sizeof(FASTX));

	if (strncmp(filename,"-",1)==0) {
		pFASTX->input = stdin;	
	} else {
		pFASTX->input = fopen(filename, "r");
		if (pFASTX->input==NULL)
			err(1, "failed to open input file '%s'", filename);
	}

	strncpy(pFASTX->input_file_name, filename, sizeof(pFASTX->input_file_name)-1);

	pFASTX->allow_input_filetype = allowed_input_filetype;
	pFASTX->allow_lowercase = allow_lowercase;
	pFASTX->allow_N = ((allow_bases & ALLOW_N)!=0) ;
	pFASTX->allow_U = ((allow_bases & ALLOW_U)!=0) ;
	pFASTX->fastq_ascii_quality_offset = fastq_ascii_quality_offset ;

	create_lookup_table(pFASTX);

	detect_input_format(pFASTX);
}

int open_output_file(const char* filename)
{
	int fd ;
	if (strncmp(filename,"-", 6)==0) {
		fd = STDOUT_FILENO;
	} else {
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );
		if (fd==-1)
			err(1, "Failed to create output file (%s)", filename);
	}
	return fd;
}

int open_output_compressor(FASTX __attribute__((unused)) *pFASTX, const char* filename)
{
	int fd;
	pid_t child_pid;
	int parent_pipe[2];
	if (pipe(parent_pipe)!=0)
		err(1,"pipe (for gzip) failed");
		
	child_pid = fork();
	if (child_pid>0) {
		/* The parent process */
		fd = parent_pipe[1];
		close(parent_pipe[0]);
		return fd;
	}

	/* The child process */

	//the compressor's STDIN is the pipe from the parent
	dup2(parent_pipe[0], STDIN_FILENO);
	close(parent_pipe[1]);

	//the compressor's STDOUT is the output file
	//(which can be the parent's STDOUT, too)
	fd = open_output_file(filename);
	dup2(fd, STDOUT_FILENO);
	
	//Run GZIP
	execlp("gzip","gzip",(char*)NULL);

	//Should never get here...
	err(1,"execlp(gzip) failed");

	return 0; //just to please the compiler
}


void fastx_init_writer(FASTX *pFASTX,
		const char *filename,
		OUTPUT_FILE_TYPE output_type, 
		int compress_output)
{
	int fd;

	if (pFASTX==NULL)
		errx(1,"Internal error: pFASTX==NULL (%s:%d)", __FILE__,__LINE__);
	if (pFASTX->input==NULL)
		errx(1,"Internal error: pFASTX not initialized (%s:%d)", __FILE__, __LINE__);

	pFASTX->compress_output = compress_output;
	if (pFASTX->compress_output)
		fd = open_output_compressor(pFASTX, filename);
	else	
		fd = open_output_file(filename);

	pFASTX->output = fdopen(fd,"w");
	if (pFASTX->output==NULL)
		err(1,"fdopen failed");

	switch(output_type)
	{
	case OUTPUT_FASTA:
		pFASTX->write_fastq = 0 ;
		pFASTX->output_sequence_id_prefix = '>';
		break ;

	case OUTPUT_FASTQ_ASCII_QUAL:
		if (! pFASTX->read_fastq) 
			errx(1,"Can't output FASTQ when input is FASTA.");
		pFASTX->write_fastq = 1;
		pFASTX->write_fastq_ascii = 1;
		pFASTX->output_sequence_id_prefix = '@';
		break ;
	
	case OUTPUT_FASTQ_NUMERIC_QUAL:
		if (! pFASTX->read_fastq) 
			errx(1,"Can't output FASTQ when input is FASTA.");
		pFASTX->write_fastq = 1;
		pFASTX->write_fastq_ascii = 0;
		pFASTX->output_sequence_id_prefix = '@';
		break ;

	case OUTPUT_SAME_AS_INPUT:
		pFASTX->write_fastq = pFASTX->read_fastq;

		//Assume we're writing ASCII format,
		pFASTX->write_fastq_ascii = 1 ;
		//But set this flag and the real format will be determined
		//when we actually read the FASTQ record
		pFASTX->copy_input_fastq_format_to_output = 1;

		pFASTX->output_sequence_id_prefix = (pFASTX->write_fastq) ? '@' : '>';
		break;

	default:
		errx(1, __FILE__ ":%d: Unknown output_type (%d)", 
			__LINE__, output_type ) ;
	}
}
	
int fastx_read_next_record(FASTX *pFASTX)
{
	char temp_qual[MAX_SEQ_LINE_LENGTH+1];

	temp_qual[MAX_SEQ_LINE_LENGTH] = 0;

	if (pFASTX==NULL)
		errx(1,"Internal error: pFASTX==NULL (%s:%d)", __FILE__,__LINE__);

	pFASTX->input_line_number++;
	if (fgets(pFASTX->input_sequence_id_prefix, MAX_SEQ_LINE_LENGTH, pFASTX->input) == NULL)
		return 0; //assume end-of-file, if we couldn't read the first line of the foursome

	chomp(pFASTX->name);

	// quick sanity check - 
	//   FASTQ files should start with '@' in the identifier line
	//   FASTA files should start with '>' in the identifier line
	if ( pFASTX->read_fastq && (pFASTX->input_sequence_id_prefix[0] != '@' ) )
		errx(1,"Invalid input: expecting FASTQ prefix character '@' on line %lld. Is this a valid FASTQ file?\n",
				pFASTX->input_line_number) ;
	if ( !pFASTX->read_fastq && (pFASTX->input_sequence_id_prefix[0] != '>') )  {
		//Extra friendly check, warn users if they fed us a multiline FASTA file
		if ( validate_nucleotides_string ( pFASTX->allowed_nucleotides, pFASTX->input_sequence_id_prefix ) ) 
			errx(1,"Invalid input: This looks like a multi-line FASTA file.\n" \
				"Line %lld contains a nucleotides string instead of a '>' prefix.\n" \
				"FASTX-Toolkit can't handle multi-line FASTA files.\n" \
				"Please use the FASTA-Formatter tool to convert this file into a single-line FASTA.\n", 
				pFASTX->input_line_number) ;
		
		// Otherwise, assume it's just bad input file
		errx(1,"Invalid input: expecting FASTA prefix character '>' on line %lld. Is this a valid FASTA file?\n",
				pFASTX->input_line_number) ;
	}



	//for the rest of the lines, if they don't appear, it's an error
	pFASTX->input_line_number++;

	if (fgets(pFASTX->nucleotides,  MAX_SEQ_LINE_LENGTH, pFASTX->input) == NULL) 
		errx(1,"Failed to read complete record, missing 2nd line (nucleotides), on line %lld\n",
			pFASTX->input_line_number);

	chomp(pFASTX->nucleotides);

	if (!validate_nucleotides_string(pFASTX->allowed_nucleotides, pFASTX->nucleotides)) 
		errx(1,"found invalid nucleotide sequence (%s) on line %lld\n",
				pFASTX->nucleotides,pFASTX->input_line_number);
	
	if (pFASTX->read_fastq) {
		pFASTX->input_line_number++;
		if (fgets(pFASTX->input_name2_prefix,  MAX_SEQ_LINE_LENGTH, pFASTX->input) == NULL) 
			errx(1,"Failed to read complete record, missing 3rd line (name-2), on line %lld\n",
				pFASTX->input_line_number);
		
		pFASTX->input_line_number++;
		if (fgets(temp_qual, sizeof(temp_qual), pFASTX->input) == NULL)
			errx(1,"Failed to read complete record, missing 4th line (quality), on line %lld\n",
				pFASTX->input_line_number);

		chomp(pFASTX->name2);
		chomp(temp_qual);
		
		if (strlen(temp_qual) == strlen(pFASTX->nucleotides)) {
			//Assume this is an ASCII quality score line, convert it to values
			convert_ascii_quality_score_line ( temp_qual, pFASTX ) ;
			pFASTX->read_fastq_ascii = 1 ;
		} else {
			//Assume this is a numeric quality score line, convert it to values
			convert_numeric_quality_score_line ( temp_qual, pFASTX ) ;
			pFASTX->read_fastq_ascii = 0 ;
		}

		//Copy the input format to the output format flag
		if (pFASTX->copy_input_fastq_format_to_output) {
			pFASTX->write_fastq_ascii = pFASTX->read_fastq_ascii;
		}
			

	}

	pFASTX->num_input_sequences++;
	pFASTX->num_input_reads += get_reads_count(pFASTX);

	return 1;
}

static void write_ascii_qual_string(FASTX *pFASTX, int length)
{
	int i;
	int rc;

	for (i=0; i<length; i++) {
		rc = fprintf(pFASTX->output, "%c", pFASTX->quality[i] + pFASTX->fastq_ascii_quality_offset ) ;
		if (rc<=0)
			err(1,"writing quality scores failed");
	}
	rc = fprintf(pFASTX->output, "\n");
	if (rc<=0)
		err(1,"writing quality scores failed");
}

static void write_numeric_qual_string(FASTX *pFASTX, int length)
{
	int i;
	int rc;
	for (i=0; i<length; i++) {
		rc = fprintf(pFASTX->output, "%d", pFASTX->quality[i] ) ;
		if (rc<=0)
			err(1,"writing quality scores failed");
		if (i<length-1) {
			rc = fprintf(pFASTX->output," ");
			if (rc<=0)
				err(1,"writing quality scores failed");
		}
	}
	rc = fprintf(pFASTX->output, "\n");
	if (rc<=0)
		err(1,"writing quality scores failed");
}

void fastx_write_record(FASTX *pFASTX)
{
	int len;
	int rc;

	if (pFASTX==NULL)
		errx(1,"Internal error: pFASTX==NULL (%s:%d)", __FILE__,__LINE__);

	
	rc = fprintf(pFASTX->output, "%c%s\n", 
			pFASTX->output_sequence_id_prefix,
			pFASTX->name ) ;
	if (rc<=0)
		err(1,"writing sequence identifier failed");
	
	rc = fprintf(pFASTX->output, "%s\n", pFASTX->nucleotides);
	if (rc<=0)
		err(1,"writing nucleotides failed");

	if (pFASTX->write_fastq) {
		rc = fprintf(pFASTX->output, "+%s\n", pFASTX->name2 ) ;
		if (rc<=0)
			err(1,"writing 2nd sequence identifier failed");

		len = strlen(pFASTX->nucleotides);	
		if (pFASTX->write_fastq_ascii)
			write_ascii_qual_string(pFASTX, len);
		else
			write_numeric_qual_string(pFASTX, len);
	}

	pFASTX->num_output_sequences++;
	pFASTX->num_output_reads += get_reads_count(pFASTX);
}

int get_reads_count(const FASTX *pFASTX)
{
	char *dash = NULL ;

	//FASTQ files are never collapsed (at least not in Gordon's Galaxy)
	if (pFASTX->read_fastq)
		return 1;

	dash = strchr(pFASTX->name,'-');

	// minus character wasn't found-
	// this sequence is most probably not collapsed
	if (dash==NULL)
		return 1;

	int count = atoi(dash+1);
	if (count>0)
		return count;
	
	return 1;
}

size_t num_input_sequences(const FASTX *pFASTX)
{
	return pFASTX->num_input_sequences;
}

size_t num_input_reads(const FASTX *pFASTX)
{
	return pFASTX->num_input_reads;
}

size_t num_output_sequences(const FASTX *pFASTX)
{
	return pFASTX->num_output_sequences;
}

size_t num_output_reads(const FASTX *pFASTX)
{
	return pFASTX->num_output_reads;
}


