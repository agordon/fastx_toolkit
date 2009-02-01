#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include <algorithm>

#include "sequence_alignment.h"

using namespace std;

void  SequenceAlignmentResults::print(std::ostream& strm) const
{
	size_t delta;
	size_t index;

	strm << "Query-Alingment = " << query_alignment << endl ;
	strm << "target-Alingment= " << target_alignment << endl ;


	strm << "Score = " << score << " ("
	     << matches << " matches, "
	     << mismatches << " mismatches, "
	     << gaps << " gaps)"
	     << std::endl ;
	
	strm << "Query = " << query_sequence
	     << "(qsize " << query_size
	     << " qstart " << query_start
	     << " qend " << query_end 
	     << std::endl ;

	strm << "Target= " << target_sequence
	     << "(tsize " << target_size
	     << " tstart " << target_start
	     << " tend " << target_end 
	     << std::endl ;

	delta = max(target_start, query_start);

	if ( delta - query_start > 0 )
		strm << std::string( delta - query_start, ' ') ;
	if ( query_start > 0 )
		strm << query_sequence.substr(0, query_start) ;
	strm << query_alignment ;

	if ( query_end+1 < query_sequence.length() )
		strm << query_sequence.substr( query_end+1 ) ;
	strm << std::endl ;

	if ( delta > 0 )
		strm << std::string( delta, ' ') ;
	for (index=0; index<query_alignment.length(); index++) {
		strm << ((query_alignment[index]==target_alignment[index]) ? '|' : '*' );
	}
	strm << std::endl;

	if ( delta - target_start > 0 )
		strm <<  std::string( delta - target_start, ' ') ;
	if ( target_start > 0 )
		strm << target_sequence.substr(0, target_start);
	strm << target_alignment;

	if ( target_end+1 < target_sequence.length() )
		strm << target_sequence.substr( target_end+1 );
	strm << std::endl;

}

SequenceAlignment::SequenceAlignment ( ) :
	_gap_panelty(-5),
	_match_panelty(1),
	_mismatch_panelty(-1)

{
}


void SequenceAlignment::set_sequences(const std::string& _query, const std::string& _target)
{
	_query_sequence = _query ;
	_target_sequence = _target ;
}

void SequenceAlignment::reset_alignment_results()
{
	_alignment_results = SequenceAlignmentResults() ;
	//
	//Reset the results
	_alignment_results.query_sequence = query_sequence() ;
	_alignment_results.target_sequence = target_sequence() ;
}

const SequenceAlignmentResults& SequenceAlignment::align ( const std::string& query, const std::string& target )
{
	set_sequences ( query, target ) ;

	reset_alignment_results();

	resize_matrix ( query_sequence().length()+1, target_sequence().length()+1 ) ;

	reset_matrix( matrix_width(), matrix_height() );
	populate_matrix();
	find_optimal_alignment();

	post_process();

	return _alignment_results;
}

void SequenceAlignment::resize_matrix(size_t width, size_t height)
{
	size_t i;

	match_matrix.resize ( width );
	for (i=0;i<width;i++)
		match_matrix[i].resize(height) ;

	score_matrix.resize ( width );
	for (i=0;i<width;i++)
		score_matrix[i].resize(height) ;

	origin_matrix.resize ( width );
	for (i=0;i<width;i++)
		origin_matrix[i].resize(height) ;

	for (size_t x=0; x<width; x++)
		for(size_t y=0;y<height;y++)
			match_matrix[x][y] = 
				match_value ( query_sequence()[x-1], target_sequence()[y-1] ) ;
}

void SequenceAlignment::post_process()
{
	
}

void SequenceAlignment::print_matrix()
{
	size_t query_index ;
	size_t target_index ;

	#if 0
	printf("Match-Matrix:\n");
	printf(" - ");
	for ( target_index=1; target_index<matrix_height(); target_index++ ) 
		printf(" %c ", target_sequence()[target_index-1] );
	printf("\n");

	for ( query_index=1; query_index<matrix_width(); query_index++ ) {

		printf(" %c ", query_sequence()[query_index-1]) ;

		for ( target_index=1 ; target_index<matrix_height(); target_index++ ) {
			printf(" %c ", match_matrix[query_index][target_index]);		
		}
		printf("\n");
	}
	#endif

	printf("Score-Matrix:\n");
	printf(" - ");
	for ( target_index=1; target_index<matrix_height(); target_index++ ) 
		printf(" %8c", target_sequence()[target_index-1] );
	printf("\n");

	for ( query_index=1; query_index<matrix_width(); query_index++ ) {

		printf("%8c ", query_sequence()[query_index-1]) ;

		for ( target_index=1 ; target_index<matrix_height(); target_index++ ) {
			DIRECTION origin = origin_matrix[query_index][target_index];

			char ch = '*';
			if (origin==FROM_UPPER)
				ch = '|';
			if (origin==FROM_LEFT)
				ch = '-';
			if (origin==FROM_UPPER_LEFT)
				ch = '\\' ;
		
			printf("%7.1f%c%c", score_matrix[query_index][target_index], ch,
					( query_index == _alignment_results.query_start && target_index == _alignment_results.target_start)?'*':' ');
		}
		printf("\n");
	}
}

#if 0
void LocalSequenceAlignment::reset_matrix( size_t width, size_t height ) 
{
	size_t x,y ;

	highest_scored_query_index = 0 ;
	highest_scored_target_index = 0 ;

	for (x=0; x<width; x++) 
		score_matrix[x][0] = 0 ;
	for (y=0; y<height; y++) 
		score_matrix[0][y] = 0 ;

}

void LocalSequenceAlignment::populate_matrix ( )
{
	size_t query_index ;
	size_t target_index ;

	ssize_t highest_score = 0 ;

	for ( query_index=1; query_index<matrix_width(); query_index++ ) {
		for ( target_index=1 ; target_index<matrix_height(); target_index++ ) {
			ssize_t score = alignment_score(query_index, target_index);

			//printf("score(q=%zu,t=%zu)=%zu\n", query_index, target_index, score ) ;
			score_matrix[query_index][target_index] = (score>0) ? score : 0 ;

			//NOTE
			// not sure ">=" is strictly correct SW (might be just ">")
			if ( score > highest_score ) {
				highest_scored_query_index = query_index ;
				highest_scored_target_index = target_index ;
				highest_score = score ;
			}
		}

	}
}

void LocalSequenceAlignment::find_optimal_alignment ( ) 
{
	size_t query_index = highest_scored_query_index ;
	size_t target_index = highest_scored_target_index;

	_alignment_results.query_end = query_index-1 ;
	_alignment_results.target_end= target_index-1 ;

	_alignment_results.score = score_matrix[query_index][target_index];

	_alignment_results.matches = 0 ;
	_alignment_results.mismatches = 0 ;

	while ( query_index > 0 || target_index > 0 ) {
		if ( score_matrix[query_index][target_index]==0)
			break ;

		//go "left" in the matrix
		if ( query_index>0 &&
		     score_matrix[query_index][target_index] == score_matrix[query_index-1][target_index] + gap_panelty() ) {

			_alignment_results.target_alignment += "-" ;
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;
			query_index--;
		}
		else
		//go "up-left" in the matrix
		if ( query_index>0 && target_index>0 &&
		     score_matrix[query_index][target_index] == 
		     	score_matrix[query_index-1][target_index-1] + match_score(query_index, target_index) ) {

			_alignment_results.target_alignment += target_sequence()[target_index-1];
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;

			(query_sequence()[query_index-1] == target_sequence()[target_index-1]) ?
				(++_alignment_results.matches) : (++_alignment_results.mismatches) ;

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

	_alignment_results.query_start = query_index ;
	_alignment_results.target_start= target_index ;

	_alignment_results.query_size = query_sequence().length();
	_alignment_results.target_size= target_sequence().length();

	std::reverse(_alignment_results.target_alignment.begin(), _alignment_results.target_alignment.end());
	std::reverse(_alignment_results.query_alignment.begin(), _alignment_results.query_alignment.end());
}
#endif

void HalfLocalSequenceAlignment::set_sequences(const std::string& _query, const std::string& _target)
{
	_query_sequence  = _query + std::string( _target.length(), 'N' ); 
	_target_sequence = std::string( _query.length(), 'N' ) + _target;
}

void HalfLocalSequenceAlignment::reset_matrix( size_t width, size_t height ) 
{
	size_t x,y ;

	highest_scored_query_index = 0 ;
	highest_scored_target_index = 0 ;

	for (x=0; x<width; x++) 
		score_matrix[x][0] = 
			//gap_panelty() * (ssize_t)x ;
			(query_sequence()[x-1]=='N') ? 0 : gap_panelty() * (ssize_t)x ;
			// 0 ;

	for (y=0; y<height; y++) 
		score_matrix[0][y] = 
			//(gap_panelty() * (ssize_t)y);
			//0;
			(target_sequence()[y-1]=='N') ? 0 : gap_panelty() * (ssize_t)y ;
			//
}

void HalfLocalSequenceAlignment::populate_matrix ( )
{
	size_t query_index ;
	size_t target_index ;
	DIRECTION origin = FROM_LEFT;

	float highest_score = 0 ;

	for ( query_index=1; query_index<matrix_width(); query_index++ ) {
		for ( target_index=1 ; target_index<matrix_height(); target_index++ ) {

			float up_score     = cell_score ( query_index, target_index-1 ) + gap_panelty() ;
			float left_score   = cell_score ( query_index-1, target_index ) + gap_panelty() ; 
			float upleft_score = cell_score ( query_index-1, target_index - 1 ) + 
						match_score ( query_index, target_index );

			//printf("query_index=%d, target_index=%d,  upscore=%f, left_score=%f, upleft_score=%f\n",
			//		query_index, target_index, up_score,left_score,upleft_score );

			float score = -100000000 ;

			if ( upleft_score > score ) {
				score = upleft_score ;
				origin = FROM_UPPER_LEFT;
			}
			if ( left_score > score ) {
				score = left_score ;
				origin = FROM_LEFT ;
			}
			if ( up_score > score ) {
				score = up_score ;
				origin = FROM_UPPER ;
			}
			//printf("query_index=%d, target_index=%d,  score=%f origin=%d\n",
			//		query_index, target_index, score, origin );
			
			score_matrix[query_index][target_index] = score ;
			origin_matrix[query_index][target_index] = origin ;

			//NOTE
			// not sure ">=" is strictly correct SW (might be just ">")
			if ( score > highest_score ) {
				highest_scored_query_index = query_index ;
				highest_scored_target_index = target_index ;
				highest_score = score ;
			}
		}
	}
}

void HalfLocalSequenceAlignment::find_optimal_alignment ( ) 
{
	#if 1
	size_t query_index = highest_scored_query_index ;
	size_t target_index = highest_scored_target_index;
	#else
	size_t query_index = query_sequence().find_last_not_of('N')+1;
	size_t target_index = target_sequence().find_last_not_of('N')+1;
	#endif

	_alignment_results.query_end = query_index-1 ;
	_alignment_results.target_end= target_index-1 ;

	_alignment_results.score = score_matrix[query_index][target_index];

	_alignment_results.matches = 0 ;
	_alignment_results.mismatches = 0 ;
	_alignment_results.gaps = 0 ;
	_alignment_results.score = 0 ;

	//printf ( "backtrace starting from (qindex=%d, tindex=%d, score=%f)\n",
	//		query_index, target_index, score_matrix[query_index][target_index]) ;

	while ( query_index > 0 && target_index > 0 ) {
		//if ( score_matrix[query_index][target_index]==0)
		//	break ;

		DIRECTION origin = origin_matrix[query_index][target_index];

		if ( origin == FROM_LEFT ) {
			_alignment_results.target_alignment += "-" ;
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;
			_alignment_results.gaps++;

			query_index--;
		}
		else
		if ( origin == FROM_UPPER_LEFT ) {
			_alignment_results.target_alignment += target_sequence()[target_index-1];
			_alignment_results.query_alignment += query_sequence()[query_index-1] ;

			(query_sequence()[query_index-1] == target_sequence()[target_index-1]) ?
				(++_alignment_results.matches) : (++_alignment_results.mismatches) ;

			query_index--;
			target_index--;
		}
		else
		if (origin == FROM_UPPER ) {
			_alignment_results.target_alignment += target_sequence()[target_index-1];
			_alignment_results.query_alignment += "-" ;
			_alignment_results.gaps++;
			target_index--;
		}
		else {
			printf("Invalid origin (%d) at query_index=%d, target_index=%d\n", 
					origin, query_index, target_index ) ;
		}
	}

	_alignment_results.query_start = query_index ;
	_alignment_results.target_start= target_index ;

	_alignment_results.query_size = query_sequence().length();
	_alignment_results.target_size= target_sequence().length();

	std::reverse(_alignment_results.target_alignment.begin(), _alignment_results.target_alignment.end());
	std::reverse(_alignment_results.query_alignment.begin(), _alignment_results.query_alignment.end());
}


void HalfLocalSequenceAlignment::post_process()
{
	//Removes the Ns which were added in 'set_sequences'
	//And adjust the results values accordingly

	//return ;
	_query_sequence.erase ( _query_sequence.find_last_not_of('N') + 1) ;
	_target_sequence.erase ( 0, _target_sequence.find_first_not_of('N') ) ;
	_alignment_results.query_sequence = _query_sequence ;
	_alignment_results.target_sequence= _target_sequence ;
	_alignment_results.query_size = _query_sequence.length();
	_alignment_results.target_size = _target_sequence.length();



	int query_n_position = _alignment_results.query_alignment.find_last_not_of('N')+1 ;
	int query_n_count    = _alignment_results.query_alignment.length() - query_n_position ;
	int target_n_count = _alignment_results.target_alignment.find_first_not_of('N') ;

	_alignment_results.query_alignment.erase( query_n_position ) ;
	_alignment_results.target_alignment.erase( 0,target_n_count ) ;

	//Update Results strucure
	_alignment_results.query_start+= target_n_count ;
	_alignment_results.query_end  -= query_n_count ;

	_alignment_results.target_start = 0;
	_alignment_results.target_end  =  _alignment_results.query_end - _alignment_results.query_start ;

	_alignment_results.query_alignment.erase ( 0, _alignment_results.query_start ) ;
	_alignment_results.target_alignment.erase ( _alignment_results.target_end+1 ) ;

	//Update match/mismatch/gap counts
	_alignment_results.matches = 0 ;
	_alignment_results.mismatches = 0 ;
	_alignment_results.gaps = 0 ;

	for (size_t index=0; index<_alignment_results.query_alignment.length(); index++) {
		char q = _alignment_results.query_alignment[index];
		char t = _alignment_results.target_alignment[index];

		if ( q == '-' || t=='-' ) {
			_alignment_results.gaps ++ ;
		} else {
			if ( q== t )
				_alignment_results.matches++;
			else
				_alignment_results.mismatches++;
		}
	}
}

