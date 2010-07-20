#ifndef __STRINGS_BUFFER_H__
#define __STRINGS_BUFFER_H__

#include <vector>
#include <tr1/memory>

class StringsBuffer
{
private:
	typedef std::tr1::shared_ptr< std::vector<char> > SHARED_STRING_BUFFER;
	SHARED_STRING_BUFFER current_buffer;
	size_t current_buffer_offset;
	size_t previous_buffer_offset;
	std::vector < SHARED_STRING_BUFFER > string_buffers;

	size_t buffer_size;

	void add_new_buffer()
	{
		SHARED_STRING_BUFFER buffer ( new std::vector<char> ) ;
		current_buffer = buffer ;
		string_buffers.push_back(buffer);

		current_buffer->resize(buffer_size);
		current_buffer_offset = 0 ;
		previous_buffer_offset = 0 ;
	}

public:
	StringsBuffer(size_t _buffer_size=128*1024*1024) :
		current_buffer_offset(0),
		previous_buffer_offset(0),
		buffer_size(_buffer_size)
	{
		add_new_buffer();
	}

	virtual ~StringsBuffer()
	{
	}

	void rollback()
	{
		current_buffer_offset = previous_buffer_offset ;
	}

	char* add_string_to_buffer(const char* str)
	{
		size_t length = strlen(str);
		if (length==0)
			return NULL;

		char* offset = (char*)allocate_buffer(length+1);
		strncpy(offset, str, length);
		offset[length]=0;
		return offset;
	}

	void* allocate_buffer0(size_t length)
	{
		void* offset = allocate_buffer(length);
		if (offset==NULL)
			return NULL;
		memset(offset, 0, length);
		return offset;

	}

	void* allocate_buffer(size_t length)
	{
		if (length==0)
			return NULL;

		if (current_buffer_offset + length > buffer_size )
			add_new_buffer();

		void* offset = (void*)&current_buffer->at(current_buffer_offset);
		previous_buffer_offset = current_buffer_offset;
		current_buffer_offset += length;
		return offset;
	}
};

#endif
