#ifndef __SEQUENCE_LIST_CONTAINER_H__
#define __SEQUENCE_LIST_CONTAINER_H__

#include <tr1/unordered_set>
#include <tr1/memory>
#include <list>
#include <vector>
#include <tr1/functional>
#include "sequence_list_loader.h"

#include <google/sparse_hash_set>

enum DUPLICATED_IDS_TYPE {
	IGNORE_DUPLICATED_IDS,
	FORBID_DUPLICATED_IDS
};

enum PAIRED_END_MARKER_TYPE {
	KEEP_PAIRED_END_MARKER,
	CHOMP_PAIRED_END_MARKER
} ;

struct eqstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	}
};



class HashedSequenceListContainer : public ISequenceIDContainer
{
	bool _allow_duplicates;
	bool _chomp_paired_end_marker;


	//typedef std::tr1::unordered_set<std::string> CONTAINER;
	typedef google::sparse_hash_set<const char*, std::tr1::hash<const char*>, eqstr> CONTAINER;

	CONTAINER container;

public:
	HashedSequenceListContainer( DUPLICATED_IDS_TYPE dup_ids, PAIRED_END_MARKER_TYPE paired_end_marker );

	virtual void add_sequence_id ( const std::string& id ) ;

	CONTAINER::const_iterator begin() const { return container.begin(); }
	CONTAINER::const_iterator end() const { return container.end(); }

	CONTAINER::iterator begin() { return container.begin(); }
	CONTAINER::iterator end() { return container.end(); }

private:
	typedef std::tr1::shared_ptr< std::vector<char> > SHARED_STRING_BUFFER;
	SHARED_STRING_BUFFER current_buffer;
	size_t current_buffer_offset;
	std::list < SHARED_STRING_BUFFER > string_buffers;

	char* add_id_to_buffer(const std::string& id);
	void add_new_id_buffer();
};

#endif
