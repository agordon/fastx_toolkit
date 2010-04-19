#ifndef __TILE_DETECTION_H__
#define __TILE_DETECTION_H__

#include <vector>
#include <list>
#include <algorithm>
#include <ostream>
#include <iomanip>

#include "adapter_hash.h"

/*
   Detection of individual tiles (k-mers)
   in the query sequence
*/
struct DetectedTile
{
	size_t query_offset ;
	size_t adapter_offset ;
	//TOOD: tile_seq is only needed for debugging and pretty-printing of information
	//      disable it with a compile-time flag?
	std::string tile_seq;
	size_t tile_size;

	DetectedTile () :
		query_offset(UINT_MAX),
		adapter_offset(UINT_MAX),
		tile_size(UINT_MAX)
	{
	}
	DetectedTile (size_t _query_offset,
			size_t _adapter_offset, const std::string &_tile_seq) :
		query_offset(_query_offset),
		adapter_offset(_adapter_offset),
		tile_seq(_tile_seq),
		tile_size(_tile_seq.length())
	{
	}

	bool acceptable_distance(const DetectedTile &other, size_t max_gap) const
	{
		int query_distance   = query_offset - other.query_offset;
		int adapter_distance = adapter_offset - other.adapter_offset ;

		return ((size_t)abs(query_distance-adapter_distance)<=max_gap);
	}
};
typedef std::vector<DetectedTile> detect_tiles_vector_type;

bool detect_tiles ( const std::string& query, SequenceHash adapter, detect_tiles_vector_type & /*output*/ detected_tiles );

inline
std::ostream& operator<< ( std::ostream& strm, const DetectedTile& tile )
{
	strm << " tile " << tile.tile_seq
		<< " query_offset= " << std::setw(3) << tile.query_offset
		<< " adapter_offset= " << std::setw(3) << tile.adapter_offset
		<< std::endl;
	return strm;
}

/*
   Consolidation of tiles into matching regions
 */
struct DetectedTileRegion
{
	size_t query_start;
	size_t query_end;

	size_t adapter_start;
	size_t adapter_end;

	size_t num_consecutive_tiles;
	size_t num_bases_covered;

	DetectedTileRegion ( ):
		query_start(UINT_MAX),
		query_end(UINT_MAX),
		adapter_start(UINT_MAX),
		adapter_end(UINT_MAX),
		num_consecutive_tiles(UINT_MAX),
		num_bases_covered(UINT_MAX)
	{
	}

	DetectedTileRegion ( size_t q_start, size_t q_end,
				size_t a_start, size_t a_end,
				size_t consecutive_tiles,
				size_t bases_covered ):
		query_start(q_start),
		query_end(q_end),
		adapter_start(a_start),
		adapter_end(a_end),
		num_consecutive_tiles(consecutive_tiles),
		num_bases_covered(bases_covered)
	{
	}

	void set_start_positions(const DetectedTile& tile)
	{
		query_start = tile.query_offset ;
		adapter_start = tile.adapter_offset;
		query_end = UINT_MAX;
		adapter_end = UINT_MAX;
		num_consecutive_tiles = 1 ;
		num_bases_covered = UINT_MAX;
	}

	void add_end_positions(const DetectedTile& tile)
	{
		query_end = tile.query_offset + tile.tile_size ;
		adapter_end = tile.adapter_offset + tile.tile_size ;
		num_consecutive_tiles++;
	}

	size_t num_tiles() const { return num_consecutive_tiles; }

	bool acceptable_distance(const DetectedTileRegion &other, size_t max_gap) const
	{
		int query_distance   = other.query_start - query_end ;
		int adapter_distance = other.adapter_start - adapter_end ;

		/*
		std::cerr << "query_distace = " << query_distance
			<< "adapter_distance = " << adapter_distance
			<< std::endl;
			*/

		return ((size_t)abs(query_distance-adapter_distance)<=max_gap);
	}

	void merge_with(const DetectedTileRegion& other)
	{
		query_start = std::min(query_start, other.query_start);
		query_end   = std::max(query_end, other.query_end);

		adapter_start = std::min(adapter_start,other.adapter_start);
		adapter_end = std::max(adapter_end,other.adapter_end);

		num_consecutive_tiles += other.num_consecutive_tiles;
	}

};

inline
std::ostream& operator<< ( std::ostream& strm, const DetectedTileRegion& tile )
{
	strm
		<< " query-start " << std::setw(3) << tile.query_start
		<< " query-end " << std::setw(3) << tile.query_end
		<< " adapter-start " << std::setw(3) << tile.adapter_start
		<< " adapter-end " << std::setw(3) << tile.adapter_end
		<< " consec-tiles " << std::setw(3) << tile.num_consecutive_tiles
/*		<< " based-covered " << tile.num_bases_covered*/
		<< std::endl;
	return strm;
};

typedef std::list<DetectedTileRegion> detected_tile_region_vector_type;

bool detect_tile_regions ( const detect_tiles_vector_type &tiles,
			detected_tile_region_vector_type & /*output*/ regions );

bool join_close_regions ( detected_tile_region_vector_type & /* in / output*/ regions );

DetectedTileRegion get_longest_region ( const detected_tile_region_vector_type & regions ) ;

#endif
