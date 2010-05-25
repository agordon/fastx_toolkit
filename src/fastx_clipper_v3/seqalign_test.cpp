#include <cstdlib>
#include <err.h>
#include <string>
#include <iostream>
#include <getopt.h>
#include "sequence_alignment.h"

using namespace std;

bool print_aligned_sequences = false;
bool print_alignment_matrix = false;
bool print_aligned_coordinates = false;
bool print_score = false;
bool local_alignment = true;
int repeat_loop = 1;
string seq1;
string seq2;

void die(const string& error)
{
	cerr << error << endl;
	exit(1);
}

void usage()
{
	cerr
<< "usage: compare_two_sequences [-a] [-m] [-g] [-l] SEQ1 SEQ2" << endl
<< "" << endl
<< "  -g  = global-alignment" << endl
<< "  -l  = local-alignment (the default)" << endl
<< "  -a  = print aligned sequences" << endl
<< "  -m  = print alignment matrix" << endl
;
	exit(1);
}

void parse_command_line(int argc, char *argv[])
{
	int c;
	while ( (c=getopt(argc, argv, "t:caglsmh")) != -1 ) {
		switch(c)
		{
		case 'a':
			print_aligned_sequences = true;
			break;
		case 'm':
			print_alignment_matrix = true;
			break;
		case 'g':
			local_alignment = false;
			break;
		case 'l':
			local_alignment = true;
			break;
		case 'c':
			print_aligned_coordinates = true;
			break;
		case 's':
			print_score = true;
			break;
		case 't':
			repeat_loop = atoi(optarg);
			if (repeat_loop<=1)
				errx(1,"'-t' argument requires a number>1");
			break;

		case 'h':
			usage();
			break;

		default:
			errx(1,"Unknown command line argument '%c'", optopt);
			break;
		}

	}
	if (optind+1>=argc)
		errx(1,"Missing two sequence arguments. See usage with '-h'");
	seq1 = argv[optind];
	seq2 = argv[optind+1];
}

int main(int argc, char* argv[])
{
	parse_command_line(argc,argv);

	SequenceAlignment sa;
	sa.set_alignment_type(local_alignment);

	for (int i=0; i<repeat_loop; ++i)
		sa.align(seq1,seq2);

	if (print_alignment_matrix)
		sa.debug_print_matrix(cout);

	if (print_score)
		sa.alignment_results().print_score(cout);
	if (print_aligned_coordinates)
		sa.alignment_results().print_aligned_coordinates(cout);
	if (print_aligned_sequences)
		sa.alignment_results().print_alignment(cout);

	return 0;
}

