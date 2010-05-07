#ifndef __LOCAL_ALIGNMENT_H__
#define __LOCAL_ALIGNMENT_H__

#include <vector>
#include <string>
#include <ostream>

struct SequenceAlignmentResults
{
	/* alignment results */
	std::string query;
	std::string target;
	std::string query_aligned;
	std::string target_aligned;

	size_t match_count;
	size_t mismatch_count;
	size_t gap_count;
	ssize_t score;

	size_t query_start;
	size_t query_end;
	size_t target_start;
	size_t target_end;

	SequenceAlignmentResults() :
		match_count(0),
		mismatch_count(0),
		gap_count(0),
		score(0)
	{
	}

	std::string aligned_mismatch_string() const;

	void print_alignment(std::ostream& strm) const;
	void print_aligned_coordinates(std::ostream &strm) const;
	void print_score(std::ostream &strm) const;
};

class SequenceAlignment
{
	ssize_t match_panelty ;
	ssize_t mismatch_panelty ;
	ssize_t gap_open_panelty ;
	/*ssize_t gap_extend_panelty ;*/

	int highest_score;
	int highest_score_i;
	int highest_score_j;
	bool smith_waterman;

	std::vector< std::vector<ssize_t> > matrix;
	std::string seq1;
	std::string seq2;

	SequenceAlignmentResults results;

public:
	SequenceAlignment();
	SequenceAlignment(ssize_t match_score,
			  ssize_t mismatch_score,
			  ssize_t gap_open_panely
			  /*ssize_t gap_extend_panelty*/ );

	void set_alignment_type(bool local_alignment) { smith_waterman = local_alignment; }

	void align(const std::string& query, const std::string &target);

	const SequenceAlignmentResults& alignment_results() const { return results; }

	void debug_print_matrix(std::ostream &strm) const;
private:
	void build_matrix();
	void find_optimal_alignment(SequenceAlignmentResults /*output*/ &results) const;

	int alignment_score(int i, int j) ;
	int match_score(int i, int j) const
	{
		return (seq1[i]==seq2[j])?match_panelty:mismatch_panelty;
	}

};

#endif
