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
#include "chomp.h"

/*
	Chomp - 
		Removes CR/LF from given string.
	
	Input - 
		string - NULL terminated string.
			 WILL BE MODIFIED!
	Output - 
		None
		
	Remarks - 
		The first CR (ASCII 13) or LF (ASCII 10) found in the string will be replaced with a NULL - 
		Effectively chomping the string.
*/
void chomp(char *string)
{
	while (*string != 0) {
		if (*string==13 || *string==10) {
			*string = 0 ;
			return;
		}
		string++;
	}
	return ;
}


