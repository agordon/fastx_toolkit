/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2010  A. Gordon (gordon@cshl.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <error.h>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <istream>
#include <getopt.h>
#include <tr1/memory>
#include <tr1/unordered_set>
#include <cstring>
#include <err.h>
#include <errno.h>

#include <gtextutils/generic_input_stream.h>

#include <file_type_detector.h>
#include <sequence_list_loader.h>
#include <sequence_list_container.h>

using namespace std;
using namespace std::tr1;

string input_filename;
string list_filename;
bool verbose = false;
size_t tabular_list_id_column=0;
bool debug_print_ids=false;

LISTTYPE input_list_file_type = TYPE_AUTO_DETECT;

struct option filter_options[] = {
	{"list",	1,	NULL,	'L'},
	{"column",	1,	NULL,	'C'},
	{"listtype",	1,	NULL,	'T'},
	{"debug",	1,	NULL,	'D'},
	{"help",	0,	NULL,	'h'},
	{"verbose",	0,	NULL,	'v'},
	{NULL,0,0,0},
};

void show_help()
{
	cout << "" << endl;
	exit(0);
}

void parse_command_line(int argc, char* argv[])
{
	int c;
	int option_index;

	while ( (c=getopt_long(argc,argv,"D:i:L:T:C:hv", filter_options, &option_index)) != -1 ) {
		switch(c)
		{
			/* Standard Options */
		case 'i':
			input_filename = optarg;
			break;

		case 'h':
			show_help();
			break;

		case 'v':
			verbose = true ;
			break ;

		case 'L':
			list_filename = optarg;
			break;

		case 'T':
			if (strlen(optarg)==0)
				errx(1,"Parameter error: missing list type (sam/bam/fastx/tabular/txt/auto)");
			if (strcasecmp(optarg,"bam")==0)
				input_list_file_type = TYPE_BAM;
			else if (strcasecmp(optarg,"sam")==0)
				input_list_file_type = TYPE_SAM;
			else if (strncasecmp(optarg,"tab",3)==0)
				input_list_file_type = TYPE_TABULAR;
			else if (strncasecmp(optarg,"fast",4)==0)
				input_list_file_type = TYPE_FASTX;
			else if (strcasecmp(optarg,"txt")==0)
				input_list_file_type = TYPE_SIMPLE_TEXT;
			else if (strcasecmp(optarg,"auto")==0)
				input_list_file_type = TYPE_AUTO_DETECT;
			else
				errx(1,"Parameter error: invalid list type '%s'", optarg);
			break;

		case 'C':
			if (strlen(optarg)==0)
				errx(1,"Parameter error: column value (for tabular list file) requires a valid number>0");
			tabular_list_id_column = atoi(optarg);

			if (tabular_list_id_column<=0)
				errx(1,"Parameter error: Invalid column value (for tabular list file) '%s'", optarg);
			break;

		case 'D':
			if (strlen(optarg)==0)
				errx(1,"Parameter error: --debug requires a valid debug option");
			if (strcasecmp(optarg,"print_ids")==0)
				debug_print_ids = true;
			else
				errx(1,"Parameter error: unknown debug option '%s'", optarg);
			break;

		default:
			exit(1);
			break;
		}
	}

	if (list_filename.empty())
		errx(1,"Parameters error: missing list filename (-L or --list) parameter");
	if (!file_type_is_readable(list_filename))
		error(1,errno,"Error: failed to open list file (%s)", list_filename.c_str());

	//if no file name specified with "-i" and there's an extra argument - assume it is the file name
	if ( input_filename.empty() && optind < argc ) {
		input_filename = argv[optind];
	}
}

int main(int argc, char* argv[])
{
	ios::sync_with_stdio(false);

	parse_command_line(argc, argv);

	//Load IDs
	HashedSequenceListContainer ids(IGNORE_DUPLICATED_IDS, CHOMP_PAIRED_END_MARKER);
	load_sequence_ids( input_list_file_type, list_filename, tabular_list_id_column, &ids);

	if (debug_print_ids) {
		copy(ids.begin(),ids.end(), ostream_iterator<string>(cout,"\n"));
		exit(0);
	}

	cerr << "Sorry, this program is not finished yet. Press ENTER to exit." << endl;
	string dummy;
	getline(cin, dummy);

	exit(1);

	return 0;
}

