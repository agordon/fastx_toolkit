#include <string>
#include <vector>
#include <algorithm>

#include "sequence_alignment.h"

using namespace std;

void  SequenceAlignmentResults::print() const
{
	printf("query  = %s\n", query_alignment.c_str()); 
	printf("target = %s\n", target_alignment.c_str()); 
}

SequenceAlignment::SequenceAlignment ( ) :
	_gap_panelty(-10),
	_match_panelty(1),
	_mismatch_panelty(-1)

{
}

const SequenceAlignmentResults& SequenceAlignment::align ( const std::string& query, const std::string& target )
{
	_query_sequence = query ;
	_target_sequence = target ;

	resize_matrix ( query.length()+1, target.length()+1 ) ;

	reset_matrix( matrix_width(), matrix_height() );
	populate_matrix();
	find_optimal_alignment();

	return _alignment_results;
}

void SequenceAlignment::resize_matrix(size_t width, size_t height)
{
	size_t i;
	matrix.resize ( width );
	for (i=0;i<width;i++)
		matrix[i].resize(height) ;
}

void SequenceAlignment::print_matrix()
{
	size_t query_index ;
	size_t target_index ;

	printf(" - ");
	for ( target_index=1; target_index<target_sequence().length(); target_index++ ) 
		printf(" %c ", target_sequence()[target_index-1] );
	printf("\n");

	for ( query_index=1; query_index<query_sequence().length(); query_index++ ) {

		printf(" %c ", query_sequence()[query_index-1]) ;

		for ( target_index=1 ; target_index<target_sequence().length(); target_index++ ) {
			printf("%2d ", matrix[query_index][target_index]);		
		}
		printf("\n");
	}
}

void LocalSequenceAlignment::reset_matrix( size_t width, size_t height ) 
{
	size_t x,y ;

	highest_scored_query_index = 0 ;
	highest_scored_target_index = 0 ;

	for (x=0; x<width; x++) 
		matrix[x][0] = 0 ;
	for (y=0; y<height; y++) 
		matrix[0][y] = 0 ;
}

void LocalSequenceAlignment::LocalSequenceAlignment::populate_matrix ( )
{
	size_t query_index ;
	size_t target_index ;

	ssize_t highest_score = 0 ;

	for ( query_index=1; query_index<query_sequence().length(); query_index++ ) {
		for ( target_index=1 ; target_index<target_sequence().length(); target_index++ ) {
			ssize_t score = alignment_score(query_index, target_index);

			//printf("score(q=%zu,t=%zu)=%zu\n", query_index, target_index, score ) ;
			matrix[query_index][target_index] = (score>0) ? score : 0 ;

			//NOTE
			// not sure ">=" is strictly correct SW (might be just ">")
			if ( score >= highest_score ) {
				highest_scored_query_index = query_index ;
				highest_scored_target_index = target_index ;
				highest_score = score ;
			}
		}
	}
}

void LocalSequenceAlignment::LocalSequenceAlignment::find_optimal_alignment ( ) 
{
	size_t query_index = highest_scored_query_index ;
	size_t target_index = highest_scored_target_index;

	while ( query_index > 0 || target_index > 0 ) {
		if ( matrix[query_index][target_index]==0)
			break ;

		//go "left" in the matrix
		if ( query_index>0 &&
		     matrix[query_index][target_index] == matrix[query_index-1][target_index] + gap_panelty() ) {

			_alignment_results.target_alignment += "-" ;
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;
			query_index--;
		}
		else
		//go "up-left" in the matrix
		if ( query_index>0 && target_index>0 &&
		     matrix[query_index][target_index] == 
		     	matrix[query_index-1][target_index-1] + match_score(query_index, target_index) ) {

			_alignment_results.target_alignment += target_sequence()[target_index-1];
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;
			query_index--;
			target_index--;
		}
		else 
		//go "up" in the matrix
		{
			_alignment_results.target_alignment += target_sequence()[target_index-1];
			_alignment_results.query_alignment += "-" ;
			target_index--;
		}
	}

	std::reverse(_alignment_results.target_alignment.begin(), _alignment_results.target_alignment.end());
	std::reverse(_alignment_results.query_alignment.begin(), _alignment_results.query_alignment.end());
}

