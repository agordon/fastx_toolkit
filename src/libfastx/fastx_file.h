#ifndef __FASTX_FILE_H__
#define __FASTX_FILE_H__

#include <string>
#include "sequence.h"

ISequenceReader* create_fastx_reader(const std::string& filename, int ASCII_quality_offset);

ISequenceReaderPE* create_fastx_pe_reader(const std::string& filename1, const std::string& filename2, int ASCII_quality_offset);

#endif
