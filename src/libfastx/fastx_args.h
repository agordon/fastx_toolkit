#ifndef __FASTX_ARGS__
#define __FASTX_ARGS__

//One day this would all be OO :-)

const char* get_input_filename();
const char* get_output_filename();
int verbose_flag();
int compress_output_flag();
FILE* get_report_file();

typedef int (*parse_argument_func)(int optind, int optc, char* optarg)  ;

int fastx_parse_cmdline( int argc, char* argv[],
			 char* program_options,
			 parse_argument_func program_parse_arg ) ;


#endif

