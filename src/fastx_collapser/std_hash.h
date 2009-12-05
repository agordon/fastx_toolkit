/*
    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)

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
#ifndef __STD_HASH__
#define __STD_HASH__


/*
 * Centralized place to load std::hash_map
 *
 * GCC needs the following hacks...
 * Other compilers/systems might require different hacks
 */

#include <ext/hash_map>
#include <ext/hash_set>

namespace std
{
	using namespace __gnu_cxx;

	struct std_string_hash
	{                                                                                           
		size_t operator()( const std::string& x ) const                                           
		{                                                                                         
			//printf("std_string_hash: hashing '%s'\n", x.c_str());
			return hash< const char* >()( x.c_str() );                                              
		}                                                                                         
	};
	
	/*
	 * 'eqstr' and 'hash_map' usage is based on http://www.sgi.com/tech/stl/hash_map.html
	 */
	struct eqstr
	{
		bool operator()(const char* s1, const char* s2) const
		{
			return strcmp(s1, s2) == 0;
		}
	};

	typedef hash_map< const char*, int, hash< const char* >, eqstr > hash_map_charptr_to_int;

	typedef hash_map< string, int, std_string_hash > hash_map_string_to_int;

	typedef hash_set < string, std_string_hash > hash_set_string ;
}

#endif

