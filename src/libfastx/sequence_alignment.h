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
	#if 0
	typedef ssize_t score_type;
	#else
	typedef float score_type;
	#endif

	typedef enum {
		FROM_UPPER = 1,
		FROM_LEFT  = 2,
		FROM_UPPER_LEFT = 3,
		FROM_NOWHERE = 4
	} DIRECTION ;

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

	const SequenceAlignmentResults& results() const { return _alignment_results; }

	score_type cell_score ( size_t query_index, size_t target_index )  const
	{
		//printf("cell_score(q=%zu,t=%zu)=%f\n", query_index, target_index, score_matrix[query_index][target_index]) ;
		return score_matrix[query_index][target_index] ;
	}

	char match_value ( const char q, const char t ) const
	{
		if ( q=='N' || t=='N' ) 
			return 'N' ;
		
		return ( q==t ) ? 'M' : 'x' ;
	}

	score_type match_score(const size_t query_index, const size_t target_index) const
	{
		if ( query_sequence()[query_index-1]=='N' && target_sequence()[target_index-1]=='N')
			return 0.0 ;

		if ( query_sequence()[query_index-1]=='N' || target_sequence()[target_index-1]=='N')
			return neutral_panelty() ;

		return (query_sequence()[query_index-1] == target_sequence()[target_index-1]) ? 
				match_panelty() : mismatch_panelty() ;
	}

	void print_matrix(std::ostream& strm = std::cout);

	ssize_t alignment_score(const size_t query_index, const size_t target_index) const
	{
		ssize_t score = -100000000;

		//Score from the left-cell
		if ( query_index > 0 )
			if ((score_matrix[query_index-1][target_index] + gap_panelty()) > score)
				score = score_matrix[query_index-1][target_index] + gap_panelty();

		//Score from the upper-cell
		if ( target_index  > 0 ) 
			if ((score_matrix[query_index][target_index-1] + gap_panelty()) > score)
				score = score_matrix[query_index][target_index-1] + gap_panelty();

		//Score from the upper-left-cell
		if ( target_index>0 && query_index> 0) {
			if (score_matrix[query_index-1][target_index-1] + match_score(query_index,target_index) > score) 
				score = score_matrix[query_index-1][target_index-1] + match_score(query_index,target_index) ;
		}
		return score;

	}

	const SequenceAlignmentResults& align ( const std::string& query, const std::string& target ) ;

protected:
	void resize_matrix(size_t width, size_t height);

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
};

#endif

