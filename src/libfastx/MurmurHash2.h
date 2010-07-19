#ifndef __MURMURHASH_H__
#define __MURMURHASH_H__

/*
http://sites.google.com/site/murmurhash/
*/
#include <sys/types.h>

typedef u_int64_t uint64_t;

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed );

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed );
uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed );

#ifdef __cplusplus
#include <string>
#include <cstring>

struct MurmurHash2_char_ptr
{
	size_t operator()(const char* str) const
	{
		//TODO:
		//Add #IFDEFs for 32bit systems, to use a different version of murmur hash
		return MurmurHash64A(str, strlen(str),0);
	}
};

struct MurmurHash2_std_string
{
	size_t operator()(const std::string& str) const
	{
		//TODO:
		//Add #IFDEFs for 32bit systems, to use a different version of murmur hash
		return MurmurHash64A(str.data(),str.length(),0);
	}
};

#endif

#endif
