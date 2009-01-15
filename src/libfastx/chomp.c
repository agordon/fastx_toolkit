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


