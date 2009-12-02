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
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <errno.h>
#include <err.h>

#include <config.h>

#include "chomp.h"
#include "fastx.h"
#include "fastx_args.h"

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif


#define MAX_SEQUENCE_LENGTH (2000) //that's pretty arbitrary... should be enough for now

const char* usage=
"usage: fastx_quality_stats [-h] [-N] [-i INFILE] [-o OUTFILE]\n" \
"Part of " PACKAGE_STRING " by A. Gordon (gordon@cshl.edu)\n" \
"\n" \
"   [-h] = This helpful help screen.\n" \
"   [-i INFILE]  = FASTQ input file. default is STDIN.\n" \
"   [-o OUTFILE] = TEXT output file. default is STDOUT.\n" \
"   [-N]         = New output format (with more information per nucleotide/cycle).\n" \
"\n"\
"The *OLD* output TEXT file will have the following fields (one row per column):\n" \
"	column	= column number (1 to 36 for a 36-cycles read solexa file)\n" \
"	count   = number of bases found in this column.\n" \
"	min     = Lowest quality score value found in this column.\n" \
"	max     = Highest quality score value found in this column.\n"	\
"	sum     = Sum of quality score values for this column.\n" \
"	mean    = Mean quality score value for this column.\n" \
"	Q1	= 1st quartile quality score.\n" \
"	med	= Median quality score.\n" \
"	Q3	= 3rd quartile quality score.\n" \
"	IQR	= Inter-Quartile range (Q3-Q1).\n" \
"	lW	= 'Left-Whisker' value (for boxplotting).\n" \
"	rW	= 'Right-Whisker' value (for boxplotting).\n" \
"	A_Count	= Count of 'A' nucleotides found in this column.\n" \
"	C_Count	= Count of 'C' nucleotides found in this column.\n" \
"	G_Count	= Count of 'G' nucleotides found in this column.\n" \
"	T_Count	= Count of 'T' nucleotides found in this column.\n" \
"	N_Count = Count of 'N' nucleotides found in this column.\n" \
"	max-count = max. number of bases (in all cycles)\n" \
"\n" \
"\n" \
"The *NEW* output format:\n" \
"	cycle (previously called 'column') = cycle number\n" \
"	max-count\n" \
"	For each nucleotide in the cycle (ALL/A/C/G/T/N):\n"
"		count   = number of bases found in this column.\n" \
"		min     = Lowest quality score value found in this column.\n" \
"		max     = Highest quality score value found in this column.\n"	\
"		sum     = Sum of quality score values for this column.\n" \
"		mean    = Mean quality score value for this column.\n" \
"		Q1	= 1st quartile quality score.\n" \
"		med	= Median quality score.\n" \
"		Q3	= 3rd quartile quality score.\n" \
"		IQR	= Inter-Quartile range (Q3-Q1).\n" \
"		lW	= 'Left-Whisker' value (for boxplotting).\n" \
"		rW	= 'Right-Whisker' value (for boxplotting).\n" \
"\n"\
"\n";
;

FILE* outfile;
int new_output_format = 0 ;

/*
	Information for each column in the solexa file.
	("Column" here refers to the number of reads in the file, usually 36)
*/
typedef enum {
	ALL = 0,
	A = 1,
	C = 2,
	G = 3,
	T = 4,
	N = 5,
	NUC_INDEX_SIZE=6
} NUCLEOTIDE_INDEX ;

char* nucleotide_index_name[NUC_INDEX_SIZE] = { "ALL","A","C","G","T","N" } ;

//Lookup table to convert from ASCII character A/C/G/T to NUCLEOTIDE_INDEX.
//Initialized in init_values().
int nuc_to_index[256];

//Information for each nucleotide of each cycle
struct nucleotide_data
{
	int min;
	int max;
	int count;
	unsigned long long sum;

	//Instead of keeping a sorted array of all the quality values (which is needed to find the median value),
	//We keep the values in this array. similar to "Couting Sort" array in "Introduction to Algorithms", page 169.
	//Each time we encounter a quality value number (in the range of MIN_QUALITY_VALUE to MAX_QUALITY_VALUE),
	//we increment the count in the corresponding index of this array.
	int bases_values_count[QUALITY_VALUES_RANGE];
};

//Information for each cycle
struct cycle_data
{
	struct nucleotide_data nucleotide[NUC_INDEX_SIZE];
};

int sequences_count;
struct cycle_data cycles[MAX_SEQUENCE_LENGTH];
FASTX fastx;

void init_values()
{
	int i,j;

	bzero ( nuc_to_index, sizeof(nuc_to_index) ) ;
	nuc_to_index['A'] = A ;
	nuc_to_index['a'] = A ;
	nuc_to_index['C'] = C ;
	nuc_to_index['c'] = C ;
	nuc_to_index['G'] = G ;
	nuc_to_index['g'] = G ;
	nuc_to_index['T'] = T ;
	nuc_to_index['t'] = T ;
	nuc_to_index['N'] = N ;
	nuc_to_index['n'] = N ;
	
	sequences_count=0;
	bzero ( &cycles, sizeof(cycles) );

	for (i=0;i<MAX_SEQUENCE_LENGTH;i++) {
		for (j=0;j<NUC_INDEX_SIZE;++j) {
			cycles[i].nucleotide[j].min = 100 ;
			cycles[i].nucleotide[j].max = -100;
		}
	}		
}

void read_file()
{
	size_t index;
	int quality_value;
	int reads_count ;
	char nucleotide ;
	int nuc_index ;

	while ( fastx_read_next_record(&fastx) ) {

		if (strlen(fastx.nucleotides) >= MAX_SEQ_LINE_LENGTH)
			errx(1, "Internal error: sequence too long (on line %llu). Hard-coded max. length is %d",
					fastx.input_line_number, MAX_SEQ_LINE_LENGTH ) ;
		
		//for each base in the sequence...
		for (index=0; index<strlen(fastx.nucleotides); index++) {

			nucleotide = fastx.nucleotides[index];

			nuc_index = nuc_to_index[(int)nucleotide];

			//Update Nucleotides Counts
			reads_count = get_reads_count(&fastx); //if this is a collapsed FASTA file, each sequence can represent multiple reads

			cycles[index].nucleotide[ALL].count += reads_count; // total counts
			cycles[index].nucleotide[nuc_index].count += reads_count ; //per-nucleotide counts


			if (fastx.read_fastq) {
				quality_value = fastx.quality[index];

				//Update the quality statistics for all nucleotides
				cycles[index].nucleotide[ALL].min = MIN ( cycles[index].nucleotide[ALL].min, quality_value ) ;
				cycles[index].nucleotide[ALL].max = MAX ( cycles[index].nucleotide[ALL].max, quality_value ) ;
				cycles[index].nucleotide[ALL].sum += quality_value;
				cycles[index].nucleotide[ALL].bases_values_count[quality_value - MIN_QUALITY_VALUE ] += reads_count ;

				//Update the quality statistics for per nucleotides
				cycles[index].nucleotide[nuc_index].min = MIN ( cycles[index].nucleotide[nuc_index].min, quality_value ) ;
				cycles[index].nucleotide[nuc_index].max = MAX ( cycles[index].nucleotide[nuc_index].max, quality_value ) ;
				cycles[index].nucleotide[nuc_index].sum += quality_value;
				cycles[index].nucleotide[nuc_index].bases_values_count[quality_value - MIN_QUALITY_VALUE ] += reads_count ;
			}
		} 

		sequences_count++;

		//DEBUG
		//if ( (fileline-1) % 10000==0 ) { fprintf(stderr,"."); fflush(stderr) ; }
	}
}

int get_nth_value(int cycle, int nucleotide, int n)
{
	int pos;
	
	if (cycle<0 || cycle>MAX_SEQUENCE_LENGTH) {
		fprintf(stderr,"Internal error at get_nth_value, cycle=%d, nucleotide=%d\n", cycle, nucleotide);
		exit(1);
	}
	if (n==0) 
		return cycles[cycle].nucleotide[nucleotide].min;

	if (n<0 || n>=cycles[cycle].nucleotide[nucleotide].count) {
		fprintf(stderr,"Internal error at get_nth_value (cycle=%d, nucleotide=%d, n=%d), count_values[%d]=%d\n",
			cycle, nucleotide, n, cycle, cycles[cycle].nucleotide[nucleotide].count ) ;
		exit(1);
	}
	
	
	
	pos = 0 ;
	while (n > 0) {
		if (cycles[cycle].nucleotide[nucleotide].bases_values_count[pos] > n)
			break;
		n -= cycles[cycle].nucleotide[nucleotide].bases_values_count[pos];
		pos++;
		while (cycles[cycle].nucleotide[nucleotide].bases_values_count[pos]==0)
			pos++;
	}
	return pos + MIN_QUALITY_VALUE ;
}

void print_nucleotide_statistics_header(const char *nuc_name)
{
	const char* headers[] = {
		"count",
		"min",
		"max",
		"sum",
		"mean",
		"Q1",
		"med",
		"Q3",
		"IQR",
		"lW",
		"rW",
		NULL
	};

	const char** header = headers;
	while ( *header != NULL ) {
		fprintf(outfile,"\t%s_%s", nuc_name, *header);
		header++;
	}
}

/*
   print quality statistics (for one nucleotide) in new format
 */
void print_nucleotide_statistics(int cycle, int nuc_index)
{
	int Q1,Q3,IQR;
	int LeftWisker, RightWisker;

	Q1 = get_nth_value ( cycle, nuc_index, cycles[cycle].nucleotide[nuc_index].count / 4 );
	Q3 = get_nth_value ( cycle, nuc_index, cycles[cycle].nucleotide[nuc_index].count * 3 / 4 );
	IQR = Q3 - Q1 ;
	
	if ( (Q1 - IQR*3/2) < cycles[cycle].nucleotide[nuc_index].min )
		LeftWisker = cycles[cycle].nucleotide[nuc_index].min;
	else
		LeftWisker = (Q1 - IQR*3/2); //TODO - make sure there's an observed value at this point
	
	if ( (Q3 + IQR*3/2) > cycles[cycle].nucleotide[nuc_index].max )
		RightWisker = cycles[cycle].nucleotide[nuc_index].max;
	else
		RightWisker = (Q3 + IQR*3/2); //TODO - make sure there's an observed value at this point

	fprintf(outfile,"\t%d\t%d\t%d\t%lld\t",
		cycles[cycle].nucleotide[nuc_index].count,
		cycles[cycle].nucleotide[nuc_index].min,
		cycles[cycle].nucleotide[nuc_index].max,
		cycles[cycle].nucleotide[nuc_index].sum);
	
	
	fprintf(outfile,"%3.2f\t%d\t%d\t%d\t",
		((double)cycles[cycle].nucleotide[nuc_index].sum)/((double)cycles[cycle].nucleotide[nuc_index].count),
		Q1,
		get_nth_value ( cycle, nuc_index, cycles[cycle].nucleotide[nuc_index].count / 2 ),
		Q3);
	
	fprintf(outfile,"%d\t%d\t%d",
		IQR,
		LeftWisker,
		RightWisker
		);

}

void print_statistics()
{
	int nuc ;
	int cycle ;
	int max_count ;

	fprintf(outfile,"cycle\tmax_count");
	for ( nuc = 0 ;nuc < NUC_INDEX_SIZE ; ++nuc ) 
		print_nucleotide_statistics_header( nucleotide_index_name[nuc] ) ;
	fprintf(outfile,"\n");

	//Maximum number of bases (out of all cycles/columns).
	//it is always equal to the count of the first column
	//(since all reads have a base at the first column, 
	// but some might not have base at later columns (if they were clipped) )
	max_count = cycles[0].nucleotide[ALL].count ;

	for (cycle=0;cycle<MAX_SEQUENCE_LENGTH;++cycle) {
		if (cycles[cycle].nucleotide[ALL].count==0)
			break;

		//Cycle number and max_count
		fprintf(outfile,"%d\t%d", cycle+1, max_count );

		for ( nuc = 0 ;nuc < NUC_INDEX_SIZE ; ++nuc ) 
			print_nucleotide_statistics( cycle, nuc ) ;
		fprintf(outfile,"\n");
	}
}

/*
   print quality statistics in old format
 */
void print_old_statistics()
{
	int i;
	int Q1,Q3,IQR;
	int LeftWisker, RightWisker;
	
	//Fields:
	fprintf(outfile,"column\t");
	fprintf(outfile,"count\tmin\tmax\tsum\t");
	fprintf(outfile,"mean\tQ1\tmed\tQ3\t");
	fprintf(outfile,"IQR\tlW\trW\t");
	fprintf(outfile,"A_Count\tC_Count\tG_Count\tT_Count\tN_Count\t");
	fprintf(outfile,"Max_count\n");
	for (i=0;i<MAX_SEQUENCE_LENGTH;i++) {
		if (cycles[i].nucleotide[ALL].count==0)
			break;
		
		Q1 = get_nth_value ( i, ALL, cycles[i].nucleotide[ALL].count / 4 );
		Q3 = get_nth_value ( i, ALL, cycles[i].nucleotide[ALL].count * 3 / 4 );
		IQR = Q3 - Q1 ;
		
		if ( (Q1 - IQR*3/2) < cycles[i].nucleotide[ALL].min )
			LeftWisker = cycles[i].nucleotide[ALL].min;
		else
			LeftWisker = (Q1 - IQR*3/2); //TODO - make sure there's an observed value at this point
		
		if ( (Q3 + IQR*3/2) > cycles[i].nucleotide[ALL].max )
			RightWisker = cycles[i].nucleotide[ALL].max;
		else
			RightWisker = (Q3 + IQR*3/2); //TODO - make sure there's an observed value at this point

		//Column number
		fprintf(outfile,"%d\t", i+1);
		
		fprintf(outfile,"%d\t%d\t%d\t%lld\t",
			cycles[i].nucleotide[ALL].count,
			cycles[i].nucleotide[ALL].min,
			cycles[i].nucleotide[ALL].max,
			cycles[i].nucleotide[ALL].sum);
		
		
		fprintf(outfile,"%3.2f\t%d\t%d\t%d\t",
			((double)cycles[i].nucleotide[ALL].sum)/((double)cycles[i].nucleotide[ALL].count),
			Q1,
			get_nth_value ( i, ALL, cycles[i].nucleotide[ALL].count / 2 ),
			Q3);
		
		fprintf(outfile,"%d\t%d\t%d\t",
			IQR,
			LeftWisker,
			RightWisker
			);
			
		fprintf(outfile,"%d\t%d\t%d\t%d\t%d\t",
			cycles[i].nucleotide[A].count,
			cycles[i].nucleotide[C].count,
			cycles[i].nucleotide[G].count,
			cycles[i].nucleotide[T].count,
			cycles[i].nucleotide[N].count);


		//Maximum number of bases (out of all cycles/columns).
		//it is always equal to the count of the first column
		//(since all reads have a base at the first column, 
		// but some might not have base at later columns (if they were clipped) )
		fprintf(outfile,"%d\n",
			cycles[0].nucleotide[ALL].count ) ;
	}
}


int parse_program_args(int __attribute__((unused)) optind, int optc, char __attribute__((unused))* optarg)
{
	switch(optc) {
		case 'N':
			new_output_format = 1 ;
			break;

		default:
			errx(1, __FILE__ ":%d: Unknown argument (%c)", __LINE__, optc ) ;
	}
	return 1;
}



void parse_commandline(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "N", parse_program_args);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE,
		get_fastq_ascii_quality_offset() );

	if (strcmp( get_output_filename(), "-" ) == 0 ) {
		outfile = stdout;
	} else {
		outfile = fopen(get_output_filename(), "w+");
		if (outfile==NULL)	
			err(1,"Failed to create output file (%s)", get_output_filename());
	}
}


int main(int argc, char* argv[])
{
	parse_commandline(argc,argv);
	init_values();
	read_file();
	if ( new_output_format )
		print_statistics();
	else	
		print_old_statistics();
	return 0;
}
