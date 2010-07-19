#ifndef __MURMURHASH_H__
#define __MURMURHASH_H__

/*
http://sites.google.com/site/murmurhash/
*/
#include <sys/types.h>

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed );
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed );
uint64_t MurmurHash64B ( const void * key, int len, unsigned int seed );

#endif
