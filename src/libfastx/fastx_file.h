#ifndef __FASTX_FILE_H__
#define __FASTX_FILE_H__

#include <string>
#include "sequence.h"

ISequenceReader* create_fastx_reader(const std::string& filename, int ASCII_quality_offset);

#endif
