#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <err.h>
#include <stdio.h>

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
		strm << std::string( delta - query_start-1, ' ') ;
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
	populate_match_matrix();

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

	query_border.resize ( width ) ;
	target_border.resize ( height ) ;

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
}

void SequenceAlignment::populate_match_matrix()
{
	for (size_t x=0; x<matrix_width(); x++)
		for(size_t y=0;y<matrix_height();y++)
			match_matrix[x][y] = 
				match_value ( query_nucleotide(x), target_nucleotide(y) ) ;
}


void SequenceAlignment::post_process()
{
	
}

void SequenceAlignment::print_matrix(std::ostream &strm) const
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
	strm << setw(2) << left << "-" << setw(7) << "-" ;
	for ( query_index=0; query_index<matrix_width(); query_index++ ) 
		strm << setw(9) << left << query_nucleotide ( query_index ) ;
	strm << endl;
	strm << setw(2) << left << "-" << setw(7) << "-" ;
	for ( query_index=0; query_index<matrix_width(); query_index++ ) 
		strm << setw(9) << left << query_border[query_index] ;
	strm << endl;

	for ( target_index=0; target_index<matrix_height(); target_index++ ) {

		strm << setw(2) << left << target_nucleotide ( target_index ) ;
		strm << setw(6) << right << target_border[target_index] << setw(1) << " ";

		for ( query_index=0 ; query_index<matrix_width(); query_index++ ) {
			char ch ;
			switch (origin ( query_index, target_index ) ) 
			{
			case FROM_UPPER:      ch = '|' ;  break ;
			case FROM_LEFT:       ch = '-' ;  break ;
			case FROM_UPPER_LEFT: ch = '\\' ;  break ;
			case FROM_NOWHERE:    ch = '=' ;  break ;
			default:              ch = '*' ; break ;
			}

			strm << left ;
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
	//_query_sequence  = _query + std::string( _target.length(), 'N' ); 
	//_target_sequence = std::string( _query.length(), 'N' ) + _target;
	_query_sequence = _query ;
	_target_sequence = _target ;
}


void HalfLocalSequenceAlignment::reset_matrix( size_t width, size_t height ) 
{
	size_t x,y ;

	highest_scored_query_index = 0 ;
	highest_scored_target_index = 0 ;

	for (x=0; x<width; x++) {
		query_border[x] = 
			//gap_panelty() * (ssize_t)x ;
			//((query_sequence()[x-1]=='N') ? neutral_panelty() : gap_panelty()) * (ssize_t)x ;
			//((query_sequence()[x-1]=='N') ? 0 : gap_panelty()) * (ssize_t)x ;
			0 ;
	}

	for (y=0; y<height; y++) {
		target_border[y] = 
			( y <= 3 ) ? 0 : (gap_panelty() * (ssize_t)(y-3));
			//0;
			//((target_sequence()[y-1]=='N') ? 0 : gap_panelty()) * (ssize_t)y ;
			//((target_sequence()[y-1]=='N') ? neutral_panelty() : gap_panelty()) * (ssize_t)y ;
	}
			
}

void HalfLocalSequenceAlignment::populate_matrix ( )
{
	size_t query_index ;
	size_t target_index ;
	DIRECTION origin = FROM_LEFT;

	score_type highest_score = -1000000 ;
	highest_scored_query_index = 0 ;
	highest_scored_target_index = 0 ;

	for ( query_index=0; query_index<matrix_width(); query_index++ ) {
		for ( target_index=0 ; target_index<matrix_height(); target_index++ ) {

			//Note:
			// 'safe_score()' can accept negative value of -1 (and will return the border value)
			score_type up_score     = safe_score(query_index,  ((ssize_t)target_index)-1) + gap_panelty() ;
			score_type left_score   = safe_score(((ssize_t)query_index)-1,target_index )  + gap_panelty() ; 
			score_type upleft_score = safe_score(((ssize_t)query_index)-1,((ssize_t)target_index)-1) + 
						nucleotide_match_score(query_index, target_index);

			//On the diagonal line, best score can not come from upper cell
			//only from left or upper-left cells
			if ( target_index>3 && target_index-3 > query_index ) {
				left_score = -100000 ;
			}

			//printf("query_index=%d, target_index=%d,  upscore=%f, left_score=%f, upleft_score=%f\n",
			//		query_index, target_index, up_score,left_score,upleft_score );

			score_type score = -100000000 ;

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

bool HalfLocalSequenceAlignment::starting_point_close_to_end_of_sequences(const size_t query_index, const size_t target_index) const
{
	if ( (size_t)query_index  >= query_sequence().length() - 2  ||
	     (size_t)target_index >= target_sequence().length() - 2 ) {
		/* We've reach either the end of the Adapter
		 * (and the adapter is shorter than the query)
		 * Or the end of the query 
		 * (and the adapter covers up to the end of the query, and then continues on)
		 *
		 * So we can safely start the alignment from this point
		 */
		return true;
	}
	else {
		/* The adapter is not covering the query until the end.
		 */
		return false;
	}
}

#undef DEBUG_STARTING_POINT
void HalfLocalSequenceAlignment::find_alignment_starting_point(ssize_t &new_query_index, ssize_t &new_target_index) const
{
	 /*
	 * Force the alignment to start from the end of the query,
	 * find the best score at the end of the query
	 *
	 * Try (desperately) to find a match that starts at the end of the query or the end of the target/adapter)
	 */
	score_type max_score = score( matrix_width()-1, matrix_height()-1 ) ;
	for ( size_t q_index = 0 ; q_index < matrix_width(); q_index++ ) {
		for ( size_t t_index = matrix_height()-2 ; t_index < matrix_height(); t_index++ ) {
			if ( origin ( q_index, t_index ) > 0 && 
				safe_score ( q_index, t_index ) > max_score ) {
				max_score = safe_score ( q_index, t_index ) ;
				#ifdef DEBUG_STARTING_POINT
				printf("Found new max score = %f at %d,%d\n", max_score, q_index, t_index ) ;
				#endif
				new_target_index = t_index ;
				new_query_index = q_index ;
			}
		}
	}
	for ( size_t q_index = matrix_width()-2 ; q_index < matrix_width(); q_index++ ) {
		for ( size_t t_index = 0 ; t_index < matrix_height(); t_index++ ) {
			if ( origin ( q_index, t_index ) > 0 && 
				safe_score ( q_index, t_index ) > max_score ) {
				max_score = safe_score ( q_index, t_index ) ;
				#ifdef DEBUG_STARTING_POINT
				printf("Found new max score = %f at %d,%d\n", max_score, q_index, t_index ) ;
				#endif
				new_target_index = t_index ;
				new_query_index = q_index ;
			}
		}
	}
	#ifdef DEBUG_STARTING_POINT
	printf("Forcing alignment from query_index=%d, target_index=%d, score=%f, origin=%d\n",
		new_query_index, new_target_index,
		score ( new_query_index, new_target_index ),
		origin ( new_query_index, new_target_index ) );
	#endif
}


#undef DEBUG_FIND_OPTIMAL_ALIGNMENT
SequenceAlignmentResults HalfLocalSequenceAlignment::find_optimal_alignment_from_point ( const size_t query_start, const size_t target_start ) const
{
	SequenceAlignmentResults results;

	results.query_sequence = query_sequence();
	results.target_sequence= target_sequence();

	ssize_t query_index = query_start;
	ssize_t target_index = target_start ;

	results.query_end = query_index ;
	results.target_end= target_index ;

	#ifdef DEBUG_FIND_OPTIMAL_ALIGNMENT
	printf ( "backtrace starting from (qindex=%d, tindex=%d, score=%f)\n",
			query_index, target_index, score_matrix[query_index][target_index]) ;
	#endif
	
	while ( query_index >= 0 && target_index >= 0 ) {
		
		const char q_nuc = query_nucleotide(query_index);
		const char t_nuc = target_nucleotide(target_index);

		const DIRECTION current_origin = origin(query_index, target_index);
		const char current_match = match ( query_index, target_index ) ;


		#ifdef DEBUG_FIND_OPTIMAL_ALIGNMENT
		const score_type current_score = score(query_index, target_index);
		printf("query_index=%d   target_index=%d  query=%c target=%c score_matrix=%3.1f origin=%d  accumulated_score = %3.2f\n",
			query_index, target_index, 
			q_nuc, t_nuc,
			current_score, 
			current_origin,
			results.score) ;
		#endif
	
		results.query_start = query_index ;
		results.target_start= target_index ;

		switch ( current_origin )
		{
		case FROM_LEFT:
			results.target_alignment += "-" ;
			results.query_alignment += q_nuc ;
			results.gaps++;
			results.score += gap_panelty();

			query_index--;
			break ;

		case FROM_UPPER_LEFT:
			results.target_alignment += t_nuc;
			results.query_alignment += q_nuc ;

			switch ( current_match ) 
			{
			case 'N':
				results.neutral_matches++ ;
				results.score += neutral_panelty();
				break ;

			case 'M':
				results.matches++;
				results.score += match_panelty();
				break;

			case 'x':
				results.mismatches++;
				results.score += mismatch_panelty();
				break ;

			default:
				errx(1,"Internal error: unknown match type (%c) at query_index=%zu, target_index=%zu\n",
					current_match, query_index, target_index ) ;
			}

			query_index--;
			target_index--;
			break ;

		case FROM_UPPER:
			results.target_alignment += t_nuc ;
			results.query_alignment += "-" ;
			results.gaps++;
			results.score += gap_panelty();

			target_index--;
			break;

		case FROM_NOWHERE:
		default:
			print_matrix();
			printf("Invalid origin (%d) at query_index=%zu, target_index=%zu\n", 
					current_origin, query_index, target_index ) ;
			printf("Query = %s\n", query_sequence().c_str());
			printf("Target= %s\n", target_sequence().c_str());
			exit(1);
		}
	}

	results.query_size = query_sequence().length();
	results.target_size= target_sequence().length();

	std::reverse(results.target_alignment.begin(),results.target_alignment.end());
	std::reverse(results.query_alignment.begin(), results.query_alignment.end());

	return results;
}

void HalfLocalSequenceAlignment::find_optimal_alignment ( ) 
{
	SequenceAlignmentResults results ;
	

	//Try to find a good alignment, 
	//starting from the highest score cell.
	results = find_optimal_alignment_from_point ( highest_scored_query_index,
						      highest_scored_target_index ) ;

	//Some heuristics:
	//If the adapter matched 7 nucleotides anywhere in the query
	//without mismatches/gaps, accept it.
	if ( results.matches >= 7 
	     && 
	     results.mismatches == 0
	     &&
	     results.gaps == 0 ) {
	     	
		_alignment_results = results ;
		return ;
	}

	if ( starting_point_close_to_end_of_sequences ( highest_scored_query_index,
						        highest_scored_target_index ) ) {
		//We're already very close to the end of the target or query,
		//can't improve much else, so return what we've got.
		_alignment_results = results ;
		return ;
	}


	//More heuristics:
	/* The adapter is not covering the query until the end.
	 * Force the alignment to start from the end of the query,
	 * find the best score at the end of the query
	 *
	 * Try (desperately) to find a match that starts at the end of the query or the end of the target/adapter)
	 */
	ssize_t query_index = highest_scored_query_index ;
	ssize_t target_index = highest_scored_target_index;
	find_alignment_starting_point ( query_index, target_index ) ;

	_alignment_results = results ;
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

