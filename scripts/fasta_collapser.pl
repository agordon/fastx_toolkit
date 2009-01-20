#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use PerlIO::gzip;

if (scalar @ARGV < 2) {
	print<<END;
	
Fasta File Collapser
Reads a Fasta file, collapses equal sequences into one line
version 0.02 (May 2008), by Assaf Gordon <gordon\@cshl.edu>

usage:
	$0 [-q] [-g] INPUT.FA OUTPUT.FA
	
	INPUT.FA  - Input file (fasta format). Can be a GZIPped file.
	OUTPUT.FA - Will be created.
	[-q]      - Quiet mode. don't write summary to STDOUT (useful for pipeline scripting)
	[-g]      - Write GZIPped output
	
Example Input File: (Sequence "ATAT" appears multiple times)
	>CSHL_3_FC0042GGAA_1_1_537_759
	ATAT
	>CSHL_3_FC0042GGAA_1_1_774_520
	TGGC
	>CSHL_3_FC0042GGAA_1_1_742_502
	ATAT
	>CSHL_3_FC0042GGAA_1_1_781_514
	TGAG
	>CSHL_3_FC0042GGAA_1_1_757_487
	TTCA
	>CSHL_3_FC0042GGAA_1_1_903_769
	ATAT
	>CSHL_3_FC0042GGAA_1_1_605_414
	TGCG
	>CSHL_3_FC0042GGAA_1_1_724_499
	ATAT

Example Output file:
	>1-4
	ATAT
	>2-1
	TGGC
	>3-1
	TGAG
	>4-1
	TTCA
	>5-1
	TGCG

END
	exit 0;
}

our $opt_q=0;
our $opt_g=0;
getopts('gq');

open(IN, "<:gzip(autopop)", "$ARGV[0]") or die "Cannot open $ARGV[0] for reading.\n";
if ($opt_g==1) {
	open(OUT, ">:gzip", "$ARGV[1]") or die "Cannot open $ARGV[1] for writing (with GZIP).\n";
} else {
	open(OUT, ">$ARGV[1]") or die "Cannot open $ARGV[1] for writing.\n";
}

my $in_count = 0;
my $out_count = 0 ;
my %sequences ;
while ( my $name = <IN> ) {
	my $sequence = <IN> ;
	
	chomp $name;
	my $prefix = substr($name,0,1,'');
	
	die "Bad input, expecting prefix character '>' got \"$name\", stopped " unless $prefix eq '>';
	die "bad input data (expected C/G/T/A sequence, got \"$sequence\"), stopped" unless $sequence =~ m/^[GTCAN]+$/ ;
	
	$sequences{$sequence}++;
	$in_count++;
}

my $index = 1 ;
while ( my ($sequence, $count) = each %sequences) {

	print OUT ">${index}-${count}\n";
	print OUT $sequence ;
	$out_count++;
	$index++;
}

if ($opt_q==0) {
	print "Collapsed $in_count sequences into $out_count sequences\n";
}
