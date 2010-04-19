#include <err.h>
#include <iostream>
#include <string>
#include "tile_detection.h"

using namespace std;


bool detect_tiles ( const std::string& query, SequenceHash adapter, detect_tiles_vector_type & /*output*/ detected_tiles )
{
	int k_step = 1;
	size_t tile_size = adapter.get_kmer_size();
	size_t query_offset ;

	detected_tiles.clear();

	for (query_offset=0 ; query_offset  < query.length() - tile_size ;
			query_offset += k_step) {
		const std::string query_kmer = query.substr(query_offset,tile_size);

		if (adapter.kmer_exists(query_kmer)) {
			size_t adapter_offset = adapter.get_kmer_offset(query_kmer);
			detected_tiles.push_back (
			       DetectedTile ( query_offset, adapter_offset, query_kmer ) ) ;
		}
	}
	if ( detected_tiles.size()<2 )
		return false;

	return true;
}

bool detect_tile_regions ( const detect_tiles_vector_type & tiles,
			detected_tile_region_vector_type & /*output*/ regions )
{
	regions.clear();
	bool last_position_stored = false;

	detect_tiles_vector_type::const_iterator it = tiles.begin();
	DetectedTile last;
	DetectedTileRegion reg;

	last = *it++;
	reg.set_start_positions(last);
	while ( it != tiles.end() ) {
		const DetectedTile& curr = *it;

		if (curr.acceptable_distance(last, 1)) {
			//continue the current region
			reg.add_end_positions(curr);
			last_position_stored = false ;
		} else {
			//Store current region
			if (reg.num_tiles()>1)
				regions.push_back(reg);
			last_position_stored = true ;

			//Start a new region
			reg.set_start_positions(curr);
		}

		++it;
		last = curr ;
	}
	if (!last_position_stored) {
		regions.push_back(reg);
	}

	return regions.size()>0;
}

bool join_close_regions ( detected_tile_region_vector_type & /* in / output*/ regions )
{
	detected_tile_region_vector_type::iterator last_it,it = regions.begin();

	last_it = it++;
	while (it != regions.end()) {
		const DetectedTileRegion &last = *last_it;
		DetectedTileRegion & curr = *it ;

		if (last.acceptable_distance(curr,2)) {
			curr.merge_with(last);
			regions.erase(last_it);
		}
		last_it = it ;
		++it;
	}
	return true;
}

struct compare_num_consec_tiles
{
	bool operator() ( const DetectedTileRegion &a,
			const DetectedTileRegion &b )
	{
		return a.num_consecutive_tiles < b.num_consecutive_tiles;
	}
};

DetectedTileRegion get_longest_region ( const detected_tile_region_vector_type & regions )
{
	return *max_element ( regions.begin(), regions.end(),
			compare_num_consec_tiles() ) ;

}
