#ifndef __LOCAL_ALIGNMENT_H__
#define __LOCAL_ALIGNMENT_H__

#include <vector>
#include <string>

class SequenceAlignment
{
	ssize_t match_panelty ;
	ssize_t mismatch_panelty ;
	ssize_t gap_open_panelty ;
	ssize_t gap_extend_panelty ;

	int highest_score;
	int highest_score_i;
	int highest_score_j;
	bool smith_waterman;

	std::vector< std::vector<ssize_t> > matrix;
	std::string seq1;
	std::string seq2;
public:
	SequenceAlignment();
	SequenceAlignment(ssize_t match_score,
			  ssize_t mismatch_score,
			  ssize_t gap_open_panely,
			  ssize_t gap_extend_panelty);

	/*
	void print() const;
	void optimal_alignment(string& std::aligned_seq1, std::string& aligned_seq2) const;
	int alignment_score() const;
	int match_score(int i, int j) const;
	*/
private:
	void build_matrix();
	int alignment_score(int i, int j) ;
};

#endif
