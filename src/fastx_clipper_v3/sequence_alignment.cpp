#include <err.h>
#include <iostream>
#include <vector>
#include <string>
#include "sequence_alignment.h"

using namespace std;

SequenceAlignment::SequenceAlignment() :
	match_panelty(1),
	mismatch_panelty(-1),
	gap_open_panelty(-1),
	gap_extend_panelty(-2),

	highest_score(-1),
	highest_score_i(-1),
	highest_score_j(-1),
	smith_waterman(true)
{
}

SequenceAlignment::SequenceAlignment(ssize_t _match_score,
			  ssize_t _mismatch_score,
			  ssize_t _gap_open_panelty,
			  ssize_t _gap_extend_panelty) :
	match_panelty(_match_score),
	mismatch_panelty(_mismatch_score),
	gap_open_panelty(_gap_open_panelty),
	gap_extend_panelty(_gap_extend_panelty),

	highest_score(-1),
	highest_score_i(-1),
	highest_score_j(-1),
	smith_waterman(true)
{
}

int SequenceAlignment::match_score(int i, int j) const
{
	return (seq1[i]==seq2[j])?1:-1;
}


void SequenceAlignment::build_matrix()
{
	// Resize the vectors to the required matrix size
	matrix.resize( seq1.length() ) ;
	for (unsigned int i=0;i<seq1.length();i++) 
		matrix[i].resize ( seq2.length() );
	
		
	// Calculate the first column and first row
	for (unsigned int i=0;i<matrix.size();i++) 
		matrix[i][0] = smith_waterman ? 0 : space_panelty * i ;
	for (unsigned int j=0;j<matrix[0].size(); j++) 
		matrix[0][j] = smith_waterman ? 0 : space_panelty * j ;
	
	for (unsigned int i=1;i<matrix.size();i++) 
		for (unsigned int j=1;j<matrix[i].size(); j++) {
			int score = alignment_score(i,j);
			matrix[i][j] = (score>0 || !smith_waterman) ? score : 0;

			if (smith_waterman && score>highest_score) {
				highest_score_i = i;
				highest_score_j = j;
				highest_score = score;
			}
		}
}

int SequenceAlignment::alignment_score(int i, int j)
{
	int score = -100000000;
	
	if ( i > 0 )
		if ((matrix[i-1][j] + space_panelty) > score)
			score = matrix[i-1][j] + space_panelty;
	
	if ( j > 0 ) 
		if ((matrix[i][j-1] + space_panelty) > score)
			score = matrix[i][j-1] + space_panelty;
	
	if ( i>0 && j > 0) {
		if (matrix[i-1][j-1] + match_score(i,j) > score) 
			score = matrix[i-1][j-1] + match_score(i,j) ;
	}
	return score;
}

void SequenceAlignment::print() const
{
	printf("   ");
	for (unsigned int j=0;j<matrix[0].size(); j++) 
		printf("    %c", seq2.at(j) );
	
	printf("\n");
	
	for (unsigned int i=0;i<matrix.size();i++) {
		printf(" %c  ", seq1.at(i) );
		for (unsigned int j=0;j<matrix[i].size(); j++) {
			printf("%4d ", matrix[i][j] ) ;
		}
		printf( "\n" );
	}
}

int SequenceAlignment::alignment_score() const
{
	unsigned int i = matrix.size()-1;
	unsigned int j = matrix[0].size()-1;
	
	return matrix[i][j];
}

void SequenceAlignment::optimal_alignment(string& aligned_seq1, string& aligned_seq2) const
{
	unsigned int i = matrix.size()-1;
	unsigned int j = matrix[0].size()-1;

	if (smith_waterman) {
		printf("Highest score on i=%d, j=%d, score=%d\n", highest_score_i, highest_score_j, highest_score);
		i = highest_score_i ;
		j = highest_score_j ;
	}
	
	while ( i > 0 || j > 0 ) {
		if (matrix[i][j]==0)
			break;

		//go "left" in the matrix
		if ( i>0 && matrix[i][j]==matrix[i-1][j]+space_panelty ) {
			aligned_seq2 += "-" ;
			aligned_seq1 += seq1[i] ;
			i--;
		}
		else
		//go "up-left" in the matrix
		if ( j>0 && i>0 && matrix[i][j]==matrix[i-1][j-1]+match_score(i,j) ) {
			aligned_seq2 += seq2[j];
			aligned_seq1 += seq1[i];
			i--;
			j--;
		}
		else
		//go "up" in the matrix
		{
			aligned_seq2 += seq2[j] ;
			aligned_seq1 += "-" ;
			j--;
		}
	}
	
	reverse(aligned_seq1.begin(), aligned_seq1.end());
	reverse(aligned_seq2.begin(), aligned_seq2.end());
}



void die(const string& error)
{
	cerr << error << endl;
	exit(1);
}

void usage()
{
	printf("usage: compare_two_sequences SEQ1 SEQ2\n");
	exit(1);
}
