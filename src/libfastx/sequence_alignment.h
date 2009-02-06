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
#ifndef __SEQUENCE_ALIGNMENT_HEADER__
#define __SEQUENCE_ALIGNMENT_HEADER__

#include <err.h>

struct SequenceAlignmentResults
{
	int alignment_found ;

	size_t query_size ;
	size_t query_start ;
	size_t query_end ;

	size_t target_size ;
	size_t target_start ;
	size_t target_end ;

	size_t gaps;
	size_t neutral_matches ;
	size_t matches ;
	size_t mismatches ;

	float score ;

	std::string query_alignment ;
	std::string target_alignment ;

	std::string query_sequence ;
	std::string target_sequence ;

	SequenceAlignmentResults() :
		alignment_found(false),
		query_size(0),
		query_start(0),
		query_end(0),

		target_size(0),
		target_start(0),
		target_end(0),

		gaps(0),
		neutral_matches(0),	
		matches(0),
		mismatches(0),

		score(0)
	{
	} 

	void print( std::ostream& ostrm = std::cout ) const;

	virtual ~SequenceAlignmentResults() {}
} ;


class SequenceAlignment
{
protected:
	typedef float score_type;

	typedef enum {
		FROM_UPPER = 1,
		FROM_LEFT  = 2,
		FROM_UPPER_LEFT = 3,
		FROM_NOWHERE = 4
		//STOP_MARKER = 5 
	} DIRECTION ;

	std::vector < score_type > query_border ;
	std::vector < score_type > target_border ;

	std::vector< std::vector< score_type >  > score_matrix ;
	std::vector< std::vector< DIRECTION >  > origin_matrix ;
	std::vector< std::vector< char > > match_matrix ;

	score_type _gap_panelty ;
	score_type _match_panelty ;
	score_type _mismatch_panelty ;
	score_type _neutral_panelty ;

 
	SequenceAlignmentResults _alignment_results ;

	std::string _query_sequence;
	std::string _target_sequence;

public:
	SequenceAlignment ( ) ;
	virtual ~SequenceAlignment() {}

	size_t matrix_width() const { return  score_matrix.size(); }
	size_t matrix_height() const { return  score_matrix[0].size(); }

	score_type gap_panelty() const { return _gap_panelty ; }
	score_type match_panelty() const { return _match_panelty ; }
	score_type mismatch_panelty() const { return _mismatch_panelty ; }
	score_type neutral_panelty() const { return _neutral_panelty ; }

	const std::string& query_sequence() const { return _query_sequence; }
	const std::string& target_sequence() const { return _target_sequence; }

	char query_nucleotide(size_t query_index) const { return _query_sequence[query_index] ; }
	char target_nucleotide(size_t target_index) const { return _target_sequence[target_index] ; }

	const SequenceAlignmentResults& results() const { return _alignment_results; }

	char match_value ( const char q, const char t ) const
	{
		if ( q=='N' || t=='N' ) 
			return 'N' ;
		
		return ( q==t ) ? 'M' : 'x' ;
	}

	char match ( const size_t query_index, const size_t target_index) const 
	{
		return match_matrix[query_index][target_index];
	}
	DIRECTION origin (  const size_t query_index, const size_t target_index) const 
	{
		return origin_matrix[query_index][target_index];
	}

	score_type score ( const size_t query_index, const size_t target_index) const 
	{
		return score_matrix[query_index][target_index];
	}

	score_type safe_score ( const ssize_t query_index, const ssize_t target_index) const 
	{
		if (query_index==-1)
			return target_border[target_index];
		if (target_index==-1)
			return query_border[query_index];

		return score_matrix[query_index][target_index];
	}

	score_type nucleotide_match_score(const size_t query_index, const size_t target_index) const
	{
		char q = query_nucleotide(query_index);
		char t = target_nucleotide(target_index);

		if ( q=='N' && t=='N' )
			return 0.0 ;

		if ( q=='N' || t=='N' )
			return neutral_panelty() ;

		return ( q==t ) ? match_panelty() : mismatch_panelty() ;
	}

	void print_matrix(std::ostream& strm = std::cout) const;

	#if 0
	score_type calculate_alignment_score(const size_t query_index, const size_t target_index) const
	{
		score_type score = -100000000;

		/*
		score_type

		//Score from the left-cell
		if ( query_index > 0 )
			if ( (score(query_index-1,target_index) + gap_panelty()) > score)
				score = score_matrix[query_index-1][target_index] + gap_panelty();

		//Score from the upper-cell
		if ( target_index  > 0 ) 
			if ((score_matrix[query_index][target_index-1] + gap_panelty()) > score)
				score = score_matrix[query_index][target_index-1] + gap_panelty();

		//Score from the upper-left-cell
		if ( target_index>0 && query_index> 0) {
			if (score_matrix[query_index-1][target_index-1] + match_score(query_index,target_index) > score) 
				score = score_matrix[query_index-1][target_index-1] + match_score(query_index,target_index) ;
		}*/
		return score;

	}
	#endif

	const SequenceAlignmentResults& align ( const std::string& query, const std::string& target ) ;

protected:
	void resize_matrix(size_t width, size_t height);
	void populate_match_matrix();

	virtual void reset_alignment_results() ; 

	virtual void set_sequences ( const std::string& _query, const std::string &target ) ;
	virtual void reset_matrix( size_t width, size_t height ) = 0 ;
	virtual void populate_matrix ( ) = 0;
	virtual void find_optimal_alignment ( ) = 0 ;
	virtual void post_process() ;
} ;

#if 0
class LocalSequenceAlignment : public SequenceAlignment
{
protected:
	size_t highest_scored_query_index ;
	size_t highest_scored_target_index ;

public:
	virtual void reset_matrix( size_t width, size_t height )  ;
	virtual void populate_matrix ( ) ;
	virtual void find_optimal_alignment ( )  ;
};
#endif


class HalfLocalSequenceAlignment : public SequenceAlignment
{
protected:
	size_t highest_scored_query_index ;
	size_t highest_scored_target_index ;

public:
	virtual void set_sequences ( const std::string& _query, const std::string &target ) ;
	virtual void reset_matrix( size_t width, size_t height )  ;
	virtual void populate_matrix ( ) ;
	virtual void find_optimal_alignment ( )  ;
	virtual void post_process() ;

	bool starting_point_close_to_end_of_sequences(const size_t query_index, const size_t target_index) const;
	void find_alignment_starting_point(ssize_t &new_query_index, ssize_t &new_target_index) const;

	SequenceAlignmentResults find_optimal_alignment_from_point ( const size_t query_start, const size_t target_start ) const ;
};

#endif

