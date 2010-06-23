#include <cstdlib>
#include <cstddef>
#include <string>
#include <tr1/unordered_set>
#include <iostream>

#include "sequence_list_loader.h"
#include "sequence_list_container.h"

using namespace std;
using namespace std::tr1;

const size_t STRING_BUFFER_SIZE (128*1024*1024);

HashedSequenceListContainer::HashedSequenceListContainer( DUPLICATED_IDS_TYPE dup_ids, PAIRED_END_MARKER_TYPE paired_end_marker ) :
	_allow_duplicates(dup_ids==IGNORE_DUPLICATED_IDS),
	_chomp_paired_end_marker(paired_end_marker==CHOMP_PAIRED_END_MARKER)
{
}

void HashedSequenceListContainer::add_new_id_buffer()
{
	SHARED_STRING_BUFFER buffer ( new std::vector<char> ) ;
	current_buffer = buffer ;
	string_buffers.push_back(buffer);

	current_buffer->resize(STRING_BUFFER_SIZE);
	current_buffer_offset = 0 ;

	container.resize ( container.size() + STRING_BUFFER_SIZE/30 );
}

char* HashedSequenceListContainer::add_id_to_buffer(const std::string& id)
{
	if (string_buffers.empty())
	       add_new_id_buffer();
	if (current_buffer_offset + id.length()+1 > STRING_BUFFER_SIZE)
		add_new_id_buffer();

	char* offset = &current_buffer->at(current_buffer_offset);
	strncpy(offset, id.c_str(), id.length());
	current_buffer_offset += id.length()+1;

	return offset;
}

bool HashedSequenceListContainer::id_exists ( const std::string& pre_id ) const
{
	string id(pre_id);

	if (_chomp_paired_end_marker && id.size()>=2) {
		const size_t len = id.size();
		if ( id[len-2]=='/' && (id[len-1]=='1' || id[len-1]=='2'))
			id.erase(len-2);
	}

	return ( container.find(id.c_str()) != container.end() );
}

void HashedSequenceListContainer::add_sequence_id ( const std::string& pre_id )
{
	string id(pre_id);

	if (_chomp_paired_end_marker && id.size()>=2) {
		const size_t len = id.size();
		if ( id[len-2]=='/' && (id[len-1]=='1' || id[len-1]=='2'))
			id.erase(len-2);
	}

	if (!_allow_duplicates) {
		if ( container.find(id.c_str()) != container.end() ) {
			cerr << "Error: duplicated sequence ID found (" << pre_id << ")" << endl;
			exit(1);
		}
	}

	char* ptr = add_id_to_buffer(id);
	container.insert(ptr);
}
