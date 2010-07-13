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
#include <climits>
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
#include <cstring>
#include <err.h>
#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <gtextutils/generic_input_stream.h>
#include <gtextutils/generic_output_stream.h>
#include <gtextutils/string_tokenize.h>

#include "libfastx/sequence.h"
#include "libfastx/file_type_detector.h"
#include "libfastx/fastx_file.h"
#include "libfastx/tab_file.h"

using namespace std;
using namespace std::tr1;

string generate_random_nucleotides(size_t len)
{
	const char nucleotides[6] = "ACGTN";

	string result;
	result.resize(len);

	for (size_t i=0;i<len;++i)
		result[i] = nucleotides[rand()%5];

	return result;
}

bool is_nucleotide_string_for_loop(const std::string &line)
{
	if (line.empty())
		return false;

	for (size_t i=0;i<line.size();++i) {
		const char c = toupper(line[i]);
		if ( ! ((c=='A') || (c=='C') || (c=='G') || (c=='T') || (c=='N')) )
			return false;
	}
	return true;
}

bool is_nucleotide_string_data(const std::string &line)
{
	if (line.empty())
		return false;

	size_t count = line.length();
	const char *pC = line.data();
	while ( count ) {
		if ( ! ((*pC=='A') || (*pC=='C') || (*pC=='G') || (*pC=='T') || (*pC=='N')) )
			return false;
		--count;
		++pC;
	}
	return true;
}


int main(int argc, char* argv[])
{
	int string_count=1000000;
	int nucleotdes_length=76;

	ios::sync_with_stdio(false);

	if (argc>1) {
		string_count = atoi(argv[1]);
		if (string_count<=0) {
			cerr << "Parameter error: please sepcify a string-count number > 0" << endl;
			exit(1);
		}
	}

	vector<string> nucs;
	nucs.resize(string_count);

	for (int i=0;i<string_count;++i)
		nucs[i] = generate_random_nucleotides(nucleotdes_length);

	bool b=false;
	for (int i=0; i<string_count;++i)
		b |= is_nucleotide_string_for_loop(nucs[i]);

	for (int i=0; i<string_count;++i)
		b |= is_nucleotide_string_data(nucs[i]);

	cout << b << endl;

	return 0;
}
