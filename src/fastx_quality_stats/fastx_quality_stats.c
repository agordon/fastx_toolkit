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
#include <getopt.h>
#include <errno.h>
#include <err.h>

#include <config.h>

#include "chomp.h"
#include "fastx.h"
#include "fastx_args.h"

#define MAX_SEQUENCE_LENGTH (MAX_SEQ_LINE_LENGTH) // as of Nov. 2008,  110 Cycles is the max... change it as necessary

const char* usage=
"usage: fastq_qual_stat [-h] [-i INFILE] [-o OUTFILE]\n" \
"\n" \
"version " VERSION " (C) 2008 by Assaf Gordon (gordon@cshl.edu)\n" \
"   [-h] = This helpful help screen.\n" \
"   [-i INFILE]  = FASTQ input file. default is STDIN.\n" \
"   [-o OUTFILE] = TEXT output file. default is STDOUT.\n" \
"\n"\
"The output TEXT file will have the following fields (one row per column):\n" \
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
"\n";
;

FILE* outfile;

/*
	Information for each column in the solexa file.
	("Column" here refers to the number of reads in the file, usually 36)
*/
struct column_data
{
	int min;
	int max;
	unsigned long long sum;
	int count;
	int A_count;
	int C_count;
	int G_count;
	int T_count;
	int N_count;
	
	//Instead of keeping a sorted array of all the quality values (which is needed to find the median value),
	//We keep the values in this array. similar to "Couting Sort" array in "Introduction to Algorithms", page 169.
	//Each time we encounter a quality value number (in the range of MIN_QUALITY_VALUE to MAX_QUALITY_VALUE),
	//we increment the count in the corresponding index of this array.
	int bases_values_count[QUALITY_VALUES_RANGE];
};

int sequences_count;
struct column_data columns[MAX_SEQUENCE_LENGTH];
FASTX fastx;

void init_values()
{
	int i,j;
	
	sequences_count=0;
	
	for (i=0;i<MAX_SEQUENCE_LENGTH;i++) {
		
		columns[i].min = 100 ;
		columns[i].max = -100;
		columns[i].sum = 0;
		columns[i].count = 0;
		columns[i].A_count = 0;
		columns[i].C_count = 0;
		columns[i].G_count = 0;
		columns[i].T_count = 0;
		columns[i].N_count = 0;
		
		for (j=0;j<QUALITY_VALUES_RANGE;j++)
			columns[i].bases_values_count[j] = 0 ;
	}		
}



void read_file()
{
	int index;
	int quality_value;
	int reads_count ;

	while ( fastx_read_next_record(&fastx) ) {

		if (strlen(fastx.nucleotides) >= MAX_SEQ_LINE_LENGTH)
			errx(1, "Internal error: sequence too long (on line %llu). Hard-coded max. length is %d",
					fastx.input_line_number, MAX_SEQ_LINE_LENGTH ) ;
		
		//for each base in the sequence...
		for (index=0; index<strlen(fastx.nucleotides); index++) {

			if (fastx.read_fastq) {
				quality_value = fastx.quality[index];

				//Update the quality statistics
				if (columns[index].min > quality_value)
					columns[index].min  = quality_value;
				if (columns[index].max < quality_value)
					columns[index].max  = quality_value;
				columns[index].sum += quality_value;
				columns[index].bases_values_count[quality_value - MIN_QUALITY_VALUE ] ++ ;
			}

			//Update Nucleotides Counts
			reads_count = get_reads_count(&fastx); //if this is a collapsed FASTA file, each sequence can represent multiple reads
			columns[index].count += reads_count;
			
			//update the base counts statistics
			switch(fastx.nucleotides[index])
			{
				case 'A': columns[index].A_count+=reads_count; break;
				case 'C': columns[index].C_count+=reads_count; break;
				case 'T': columns[index].T_count+=reads_count; break;
				case 'G': columns[index].G_count+=reads_count; break;
				case 'N': columns[index].N_count+=reads_count; break;

				/* This shoudn't really happen, as 'fastx_read_next_record' should catch invalid values */
				default: errx(1, "Internal error: invalid base value (%c)!", fastx.nucleotides[index]) ;
			}
		} 

		sequences_count++;

		//DEBUG
		//if ( (fileline-1) % 10000==0 ) { fprintf(stderr,"."); fflush(stderr) ; }
	}
}

int get_nth_value(int base_index, int n)
{
	int pos;
	
	if (base_index<0 || base_index>MAX_SEQUENCE_LENGTH) {
		fprintf(stderr,"Internal error at get_nth_value, base_index=%d\n", base_index);
		exit(1);
	}
	if (n<0 || n>=columns[base_index].count) {
		fprintf(stderr,"Internal error at get_nth_value (base_index=%d, n=%d), count_values[%d]=%d\n",
			base_index, n, base_index, columns[base_index].count ) ;
		exit(1);
	}
	
	if (n==0) 
		return columns[base_index].min;
	
	
	pos = 0 ;
	while (n > 0) {
		if (columns[base_index].bases_values_count[pos] > n)
			break;
		n -= columns[base_index].bases_values_count[pos];
		pos++;
		while (columns[base_index].bases_values_count[pos]==0)
			pos++;
	}
	return pos + MIN_QUALITY_VALUE ;
}

void print_statistics()
{
	int i;
	int Q1,Q3,IQR;
	int LeftWisker, RightWisker;
	
	//Fields:
	fprintf(outfile,"column\t");
	fprintf(outfile,"count\tmin\tmax\tsum\t");
	fprintf(outfile,"mean\tQ1\tmed\tQ3\t");
	fprintf(outfile,"IQR\tlW\trW\t");
	fprintf(outfile,"A_Count\tC_Count\tG_Count\tT_Count\tN_Count\n");
	for (i=0;i<MAX_SEQUENCE_LENGTH;i++) {
		if (columns[i].count==0)
			break;
		
		Q1 = get_nth_value ( i, columns[i].count / 4 );
		Q3 = get_nth_value ( i, columns[i].count * 3 / 4 );
		IQR = Q3 - Q1 ;
		
		if ( (Q1 - IQR*3/2) < columns[i].min )
			LeftWisker = columns[i].min;
		else
			LeftWisker = (Q1 - IQR*3/2); //TODO - make sure there's an observed value at this point
		
		if ( (Q3 + IQR*3/2) > columns[i].max )
			RightWisker = columns[i].max;
		else
			RightWisker = (Q3 + IQR*3/2); //TODO - make sure there's an observed value at this point

		//Column number
		fprintf(outfile,"%d\t", i+1);
		
		fprintf(outfile,"%d\t%d\t%d\t%lld\t",
			columns[i].count,
			columns[i].min,
			columns[i].max,
			columns[i].sum);
		
		
		fprintf(outfile,"%3.2f\t%d\t%d\t%d\t",
			((double)columns[i].sum)/((double)columns[i].count),
			Q1,
			get_nth_value ( i, columns[i].count / 2 ),
			Q3);
		
		fprintf(outfile,"%d\t%d\t%d\t",
			IQR,
			LeftWisker,
			RightWisker
			);
			
		fprintf(outfile,"%d\t%d\t%d\t%d\t%d\n",
			columns[i].A_count,
			columns[i].C_count,
			columns[i].G_count,
			columns[i].T_count,
			columns[i].N_count);
	}
}

void parse_commandline(int argc, char* argv[])
{
	fastx_parse_cmdline(argc, argv, "", NULL);

	fastx_init_reader(&fastx, get_input_filename(), 
		FASTA_OR_FASTQ, ALLOW_N, REQUIRE_UPPERCASE);

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
	print_statistics();
	return 0;
}
