#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <err.h>

#include "sequence_alignment.h"

using namespace std;

void  SequenceAlignmentResults::print(std::ostream& strm) const
{
	size_t delta;
	size_t index;

	strm << "Query-Alingment = " << query_alignment << endl ;
	strm << "target-Alingment= " << target_alignment << endl ;


 	strm << (alignment_found ? "Alignment Found" : "Alignment NOT found") << endl;
	strm << "Score = " << score << " ("
	     << matches << " matches, "
	     << neutral_matches << " neutral-matches, "
	     << mismatches << " mismatches, "
	     << gaps << " gaps) "
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

        strm << endl;

	delta = max(target_start, query_start);


	//Spaces before the query string
	if ( delta - query_start > 0 )
		strm << std::string( delta - query_start, ' ') ;
	//Un-Aligned query part (prefix)
	if ( query_start > 0 )
		strm << query_sequence.substr(0, query_start-1) ;
	//Aligned query part
	strm << "(" << query_alignment << ")";
	//Un-Aligned query part (suffix)
	if ( query_end < query_sequence.length() )
		strm << query_sequence.substr( query_end+1 ) ;
	strm << std::endl ;

	//Alignment bars
	if ( delta > 0 )
		strm << std::string( delta-1, ' ') ;
	strm << "(" ;
	for (index=0; index<query_alignment.length(); index++) {
		strm << ((query_alignment[index]==target_alignment[index]) ? '*' : '|' );
	}
	strm << ")" ;
	strm << std::endl;

	//Spaces before the target string
	if ( delta - target_start > 0 )
		strm <<  std::string( delta - target_start, ' ') ;
	//Un-Aligned target part (prefix)
	if ( target_start > 0 )
		strm << target_sequence.substr(0, target_start-1);
	//Aligned target part
	strm << "(" << target_alignment << ")";

	//Un-Aligned target part (suffix)
	if ( target_end < target_sequence.length() )
		strm << target_sequence.substr( target_end+1 );
	strm << std::endl;

}

SequenceAlignment::SequenceAlignment ( ) :
	_gap_panelty(-5),
	_match_panelty(1),
	_mismatch_panelty(-1),
	_neutral_panelty(0.1)

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

	resize_matrix ( query_sequence().length(), target_sequence().length() ) ;

	reset_matrix( matrix_width(), matrix_height() );
	populate_matrix();
	find_optimal_alignment();

	post_process();

	return _alignment_results;
}

void SequenceAlignment::resize_matrix(size_t width, size_t height)
{
	size_t i;

	if ( matrix_width() >= width && matrix_height() >= height )
		return ;

	score_matrix.resize ( width );
	for (i=0;i<width;i++)
		score_matrix[i].resize(height) ;

	origin_matrix.resize ( width );
	for (i=0;i<width;i++)
		origin_matrix[i].resize(height) ;

	match_matrix.resize ( width );
	for (i=0;i<width;i++) {
		match_matrix[i].resize(height) ;
	}

	for (size_t x=0; x<width; x++)
		for(size_t y=0;y<height;y++)
			match_matrix[x][y] = 
				match_value ( query_nucleotide(x), target_nucleotide(y) ) ;
}

void SequenceAlignment::post_process()
{
	
}

void SequenceAlignment::print_matrix(std::ostream &strm)
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

	strm << "Score-Matrix:" << endl ;

	//Print Target nucleotides
	strm << setw(2) << left << " - " ;
	for ( target_index=0; target_index<matrix_height(); target_index++ ) 
		strm << setw(9) << left << target_nucleotide ( target_index ) ;
	strm << endl;

	for ( query_index=0; query_index<matrix_width(); query_index++ ) {

		strm << setw(2) << left << query_nucleotide ( query_index ) ;

		for ( target_index=0 ; target_index<matrix_height(); target_index++ ) {
			char ch ;
			switch (origin ( query_index, target_index ) ) 
			{
			case FROM_UPPER:      ch = '|' ;  break ;
			case FROM_LEFT:       ch = '-' ;  break ;
			case FROM_UPPER_LEFT: ch = '\\' ;  break ;
			case FROM_NOWHERE:    ch = '=' ;  break ;
			default:              ch = '*' ; break ;
			}

			strm << setw(1) << match(query_index,target_index);
			strm << setw(1) << ch ;
			strm << setw(7) << fixed << setprecision(1) 
			     << score(query_index,target_index) ;
		}
		strm << endl;
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

	for (x=0; x<width; x++) {
		score_matrix[x][0] = 
			//gap_panelty() * (ssize_t)x ;
			//((query_sequence()[x-1]=='N') ? neutral_panelty() : gap_panelty()) * (ssize_t)x ;
			//((query_sequence()[x-1]=='N') ? 0 : gap_panelty()) * (ssize_t)x ;
			0 ;
		origin_matrix[x][0] = FROM_LEFT ;
	}

	for (y=0; y<height; y++) {
		score_matrix[0][y] = 
			//(gap_panelty() * (ssize_t)y);
			0;
			//((target_sequence()[y-1]=='N') ? 0 : gap_panelty()) * (ssize_t)y ;
			//((target_sequence()[y-1]=='N') ? neutral_panelty() : gap_panelty()) * (ssize_t)y ;
		origin_matrix[0][y] = FROM_UPPER;
	}
			
}

void HalfLocalSequenceAlignment::populate_matrix ( )
{
	size_t query_index ;
	size_t target_index ;
	DIRECTION origin = FROM_LEFT;

	float highest_score = -1000000 ;
	highest_scored_query_index = -1 ;
	highest_scored_target_index = -1 ;

	for ( query_index=1; query_index<matrix_width(); query_index++ ) {
		for ( target_index=1 ; target_index<matrix_height(); target_index++ ) {
		//for ( target_index=1 ; target_index<=query_index; target_index++ ) {

			float up_score     = score(query_index,  target_index-1) + gap_panelty() ;
			float left_score   = score(query_index-1,target_index )  + gap_panelty() ; 
			float upleft_score = score(query_index-1,target_index-1) + 
						nucleotide_match_score(query_index, target_index);

			//On the diagonal line, best score can not come from upper cell
			//only from left or upper-left cells
			//if ( target_index == query_index )
			//	left_score = -100000 ;

			//printf("query_index=%d, target_index=%d,  upscore=%f, left_score=%f, upleft_score=%f\n",
			//		query_index, target_index, up_score,left_score,upleft_score );

			float score = -100000000 ;

			if ( upleft_score > score ) {
				score = upleft_score ;
				origin = FROM_UPPER_LEFT;
			}
			if ( up_score > score ) {
				score = up_score ;
				origin = FROM_UPPER ;
			}
			if ( left_score > score ) {
				score = left_score ;
				origin = FROM_LEFT ;
			}
			//printf("query_index=%d, target_index=%d,  score=%f origin=%d\n",
			//		query_index, target_index, score, origin );
			
			/*if (score<0) {
				score = 0 ;
				origin = FROM_NOWHERE ;
			}*/
			
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
	ssize_t query_index = highest_scored_query_index ;
	ssize_t target_index = highest_scored_target_index;
	#else
	size_t target_index = matrix_height()-1;

	size_t query_index = 0;
	ssize_t max_score = score_matrix[query_index][target_index] ;
	for ( size_t index = 0 ; index < matrix_width(); index++ ) {
		if ( score_matrix[index][target_index] > max_score ) {
			max_score = score_matrix[index][target_index];
			query_index = index ;
		}
	}
	#endif

	_alignment_results.query_end = query_index ;
	_alignment_results.target_end= target_index ;

	_alignment_results.score = score_matrix[query_index][target_index];

	_alignment_results.neutral_matches = 0 ;
	_alignment_results.matches = 0 ;
	_alignment_results.mismatches = 0 ;
	_alignment_results.gaps = 0 ;
	_alignment_results.score = 0 ;

	//printf ( "backtrace starting from (qindex=%d, tindex=%d, score=%d)\n",
	//		query_index, target_index, score_matrix[query_index][target_index]) ;
	
	score_type current_score = 0;

	while ( query_index >= 0 && target_index >= 0 ) {
		
		const char q_nuc = query_nucleotide(query_index);
		const char t_nuc = target_nucleotide(target_index);

		//Gordon's improvement over Alex's clipper
		if (t_nuc=='N')
			break ;

		const DIRECTION current_origin = origin(query_index, target_index);
		const char current_match = match ( query_index, target_index ) ;
		current_score = score(query_index, target_index);
		

		#if 0
		printf("query_index=%d   target_index=%d  query=%c target=%c score_matrix=%3.1f origin=%d  accumulated_score = %3.2f\n",
			query_index, target_index, 
			q_nuc, t_nuc,
			current_score, 
			current_origin,
			_alignment_results.score) ;
		#endif
	
		_alignment_results.query_start = query_index ;
		_alignment_results.target_start= target_index ;

		//Original Alex's clipper behaviour, remove for better(?) gordon behaviour
		//if ( current_score < match_panelty() )
		//	break ;

		/*
		if ( current_origin == FROM_LEFT || current_origin == FROM_UPPER_LEFT || current_origin == STOP_MARKER ) {
			_alignment_results.query_alignment += q_nuc ;
		} else {
			_alignment_results.query_alignment += "-" ;
		}

		if ( current_origin == FROM_UPPER || current_origin == FROM_UPPER_LEFT || current_origin = STOP_MARKER ) {
			_alignment_results.target_alignment += t_nuc;
		} else {
			_alignment_results.target_alignment += "-" ;
		}*/

		switch ( current_origin )
		{
		case FROM_LEFT:
			_alignment_results.target_alignment += "-" ;
			_alignment_results.query_alignment += q_nuc ;
			_alignment_results.gaps++;
			_alignment_results.score += gap_panelty();

			query_index--;
			break ;

		case FROM_UPPER_LEFT:
			_alignment_results.target_alignment += t_nuc;
			_alignment_results.query_alignment += q_nuc ;

			switch ( current_match ) 
			{
			case 'N':
				_alignment_results.neutral_matches++ ;
				_alignment_results.score += neutral_panelty();
				break ;

			case 'M':
				_alignment_results.matches++;
				_alignment_results.score += match_panelty();
				break;

			case 'x':
				_alignment_results.mismatches++;
				_alignment_results.score += mismatch_panelty();
				break ;

			default:
				errx(1,"Internal error: unknown match type (%c) at query_index=%d, target_index=%d\n",
					current_match, query_index, target_index ) ;
			}

			query_index--;
			target_index--;
			break ;

		case FROM_UPPER:
			//HalfLocal Alignment doesn't allow insertions in the query
			//break;

			_alignment_results.target_alignment += t_nuc ;
			_alignment_results.query_alignment += "-" ;
			_alignment_results.gaps++;
			_alignment_results.score += gap_panelty();

			target_index--;
			break;

	/*	case STOP_MARKER:
			_alignment_results.target_alignment += t_nuc;
			_alignment_results.query_alignment += q_nuc ;
			break ;*/

		case FROM_NOWHERE:

		default:
			print_matrix();
			printf("Invalid origin (%d) at query_index=%d, target_index=%d\n", 
					current_origin, query_index, target_index ) ;
			printf("Query = %s\n", query_sequence().c_str());
			printf("Target= %s\n", target_sequence().c_str());
			exit(1);
		}
	}

	/*_alignment_results.score = 
			_alignment_results.matches * match_panelty() +
			_alignment_results.mismatches * mismatch_panelty() +
			_alignment_results.gaps * gap_panelty() +
			_alignment_results.neutral_matches * neutral_panelty() ;*/
	
	if (query_index<=0) {
		_alignment_results.alignment_found = true ;
	}

	//_alignment_results.score = current_score ;


	_alignment_results.query_size = query_sequence().length();
	_alignment_results.target_size= target_sequence().length();

	std::reverse(_alignment_results.target_alignment.begin(), _alignment_results.target_alignment.end());
	std::reverse(_alignment_results.query_alignment.begin(), _alignment_results.query_alignment.end());
}

void HalfLocalSequenceAlignment::post_process()
{
#if 0
	//Removes the Ns which were added in 'set_sequences'
	//And adjust the results values accordingly

	//return ;
	//_query_sequence.erase ( _query_sequence.find_last_not_of('N') + 1) ;
	//_target_sequence.erase ( 0, _target_sequence.find_first_not_of('N') ) ;

	_alignment_results.query_sequence.erase ( _query_sequence.find_last_not_of('N') + 1) ;
	_alignment_results.target_sequence.erase ( 0, _target_sequence.find_first_not_of('N') ) ;
	_alignment_results.query_size = _alignment_results.query_sequence.length();
	_alignment_results.target_size = _alignment_results.target_sequence.length();


	size_t query_n_position = _alignment_results.query_alignment.find_last_not_of('N') ;
	int query_n_count;

	if ( query_n_position != string::npos )
		query_n_count = _alignment_results.query_alignment.length() - query_n_position ;
	else
		query_n_count = 0 ;

	int target_n_count = _alignment_results.target_alignment.find_first_not_of('N') ;

	if (query_n_position != string::npos )
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
#endif 
}

