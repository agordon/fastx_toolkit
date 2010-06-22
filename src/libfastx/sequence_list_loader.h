#ifndef __SEQUENCE_LIST_LOADER_H__
#define __SEQUENCE_LIST_LOADER_H__

#include <string>

typedef enum {
	TYPE_AUTO_DETECT=0,
	TYPE_TABULAR=1,
	TYPE_SAM=2,
	TYPE_BAM=3,
	TYPE_FASTX=4,
	TYPE_SIMPLE_TEXT=5
} LISTTYPE ;

class ISequenceIDContainer
{
public:
	virtual ~ISequenceIDContainer() { }
	virtual void add_sequence_id ( const std::string& id ) = 0 ;
};

void load_sequence_ids_fasta(const std::string& filename, ISequenceIDContainer *pContainer);
void load_sequence_ids_fastq(const std::string& filename, ISequenceIDContainer *pContainer);
void load_sequence_ids_txt(const std::string& filename, ISequenceIDContainer *pContainer);
void load_sequence_ids_tabular(const std::string& filename, size_t column, ISequenceIDContainer *pContainer);

void load_sequence_ids(LISTTYPE input_file, const std::string& filename, size_t optional_column, ISequenceIDContainer *pContainer);

#endif
