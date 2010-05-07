#include <err.h>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include "sequence_alignment.h"

using namespace std;

std::string SequenceAlignmentResults::aligned_mismatch_string() const
{
	string s(query_aligned.length(), ' ');
	for (size_t i=0;i<query_aligned.length();++i){
		if ( query_aligned[i] != '-'
			&& target_aligned[i] != '-'
			&& query_aligned[i] != target_aligned[i])
			s[i] = '|';
	}
	return s;
}

void SequenceAlignmentResults::print_alignment(std::ostream& strm) const
{
	strm << "query : " << query_aligned << endl
	     << "        " << aligned_mismatch_string() << endl
	     << "target: " << target_aligned << endl
	     << endl;
}

void SequenceAlignmentResults::print_aligned_coordinates(std::ostream &strm) const
{
	strm << "qStart = " << setw(3) << query_start
		<< "  qEnd = " << setw(3) << query_end
		<< endl;
	strm << "tStart = " << setw(3) << target_start
		<< "  tEnd = " << setw(3) << target_end
		<< endl;
}

void SequenceAlignmentResults::print_score(std::ostream &strm) const
{
	strm << "score = " << score
		<<  " (" << match_count << " match" << ((match_count!=1)?"es":"") << "  "
		<< mismatch_count << " mismatch" << ((mismatch_count!=1)?"es":"") << "  "
		<< gap_count << " gap" << ((gap_count!=1)?"s":"") << ")"
		<< endl ;
}

SequenceAlignment::SequenceAlignment() :
	match_panelty(1),
	mismatch_panelty(-1),
	gap_open_panelty(-1),
/*	gap_extend_panelty(-2),*/

	highest_score(-1),
	highest_score_i(-1),
	highest_score_j(-1),
	smith_waterman(true)
{
}

SequenceAlignment::SequenceAlignment(ssize_t _match_score,
			  ssize_t _mismatch_score,
			  ssize_t _gap_open_panelty
			  /*ssize_t _gap_extend_panelty*/ ) :
	match_panelty(_match_score),
	mismatch_panelty(_mismatch_score),
	gap_open_panelty(_gap_open_panelty),
/*	gap_extend_panelty(_gap_extend_panelty),*/

	highest_score(-1),
	highest_score_i(-1),
	highest_score_j(-1),
	smith_waterman(true)
{
}

void SequenceAlignment::align(const std::string& _seq1, const std::string &_seq2)
{
	seq1 = _seq1;
	seq2 = _seq2;

	build_matrix();

	find_optimal_alignment(results);
}


void SequenceAlignment::build_matrix()
{
	// Resize the vectors to the required matrix size
	matrix.resize( seq1.length()+1 ) ;
	for (unsigned int i=0;i<seq1.length()+1;i++)
		matrix[i].resize ( seq2.length()+1 );

	matrix[0][0] = 0 ;

	// Calculate the first column and first row
	for (unsigned int i=0;i<matrix.size();i++)
		matrix[i][0] = smith_waterman ? 0 : gap_open_panelty * i ;
	for (unsigned int j=0;j<matrix[0].size(); j++)
		matrix[0][j] = smith_waterman ? 0 : gap_open_panelty * j ;

	for (unsigned int i=1;i<matrix.size();i++)
		for (unsigned int j=1;j<matrix[i].size(); j++) {
			int score = alignment_score(i,j);

			//cout << "i=" << i << "  j=" << j << "  score= "<< score <<endl;

			matrix[i][j] = (score>0 || !smith_waterman) ? score : 0;

			if (smith_waterman && score>highest_score) {
				highest_score_i = i;
				highest_score_j = j;
				highest_score = score;
			}
		}
	//cout << "highest_score_i=" << highest_score_i << "  highest_score_j=" << highest_score_j << "  highest_scorescore= "<< highest_score <<endl;
}

int SequenceAlignment::alignment_score(int i, int j)
{
	int score = -100000 ;

	if ( i>0 )
		if ((matrix[i-1][j] + gap_open_panelty) > score) {
			//cout << "a" << endl;
			score = matrix[i-1][j] + gap_open_panelty;
		}

	if ( j>0 )
		if ((matrix[i][j-1] + gap_open_panelty) > score) {
			//cout << "b" << endl;
			score = matrix[i][j-1] + gap_open_panelty;
		}

	if ( i>0 && j>0) {
		if (matrix[i-1][j-1] + match_score(i-1,j-1) > score) {
			//cout << "c" << endl;
			score = matrix[i-1][j-1] + match_score(i-1,j-1) ;
		}
	}
	return score;
}

void SequenceAlignment::debug_print_matrix(std::ostream &strm) const
{
	//Print nucleotides
	strm << "        ";
	for (unsigned int j=1;j<matrix[0].size(); j++)
		strm << "   " << seq2.at(j-1);
	strm << endl;

	for (unsigned int i=0;i<matrix.size();i++) {

		//Print Numcleotide
		if (i>0)
			strm << " " << seq1.at(i-1) << "  " ;
		else
			strm << "    " ;
		for (unsigned int j=0;j<matrix[i].size(); j++) {
			strm << setw(4) << matrix[i][j];
		}
		strm << endl;
	}
}

/*
int SequenceAlignment::alignment_score() const
{
	unsigned int i = matrix.size()-1;
	unsigned int j = matrix[0].size()-1;

	return matrix[i][j];
}
*/

void SequenceAlignment::find_optimal_alignment( SequenceAlignmentResults /* output */ &results ) const
{
	unsigned int i = matrix.size()-1;
	unsigned int j = matrix[0].size()-1;

	results.match_count = 0 ;
	results.gap_count = 0 ;
	results.mismatch_count = 0 ;
	results.score = 0 ;
	results.query_aligned = "" ;
	results.target_aligned = "" ;

	if (smith_waterman) {
		//printf("Highest score on i=%d, j=%d, score=%d\n", highest_score_i, highest_score_j, highest_score);
		i = highest_score_i ;
		j = highest_score_j ;
	}

	results.query_end = i-1;
	results.target_end = j-1;

	results.score = matrix[i][j] ;

	while ( i > 0 || j > 0 ) {
		if (smith_waterman && matrix[i][j]==0)
			break;

		//go "left" in the matrix
		if ( i>0 && matrix[i][j]==matrix[i-1][j]+gap_open_panelty ) {
			results.query_aligned += seq1[i-1] ;
			results.target_aligned += "-" ;
			results.gap_count++;
			i--;
		}
		else
		//go "up-left" in the matrix
		if ( j>0 && i>0 && matrix[i][j]==matrix[i-1][j-1]+match_score(i-1,j-1) ) {
			results.query_aligned += seq1[i-1];
			results.target_aligned += seq2[j-1];
			i--;
			j--;
			if (match_score(i,j)==match_panelty)
				results.match_count++;
			else
				results.mismatch_count++;
		}
		else
		//go "up" in the matrix
		{
			results.query_aligned += "-" ;
			results.target_aligned += seq2[j-1] ;
			j--;
			results.gap_count++;
		}
	}

	results.query_start = i;
	results.target_start = j;

	std::reverse(results.query_aligned.begin(), results.query_aligned.end());
	std::reverse(results.target_aligned.begin(), results.target_aligned.end());
}


