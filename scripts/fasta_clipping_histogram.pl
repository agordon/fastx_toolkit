#!/usr/bin/perl
use strict;
use warnings;
use GD::Graph::bars;
use Data::Dumper;
use PerlIO::gzip;

if (scalar @ARGV==0) {
	print<<END;
	
Create a Linker Clipping Information Histogram

usage: $0 INPUT_FILE.FA OUTPUT_FILE.PNG

	INPUT_FILE.FA   = input file (in FASTA format, can be GZIPped)
	OUTPUT_FILE.PNG = histogram image

END
	exit 0;
}

#
# Read parameters
#
open(IN, "<:gzip(autopop)", "$ARGV[0]") or die "Cannot open input file $ARGV[0]\n";
open(OUT, ">$ARGV[1]") or die "Cannot create output file $ARGV[1]\n";
binmode OUT;

my %histogram ;

while (my $name = <IN>) {
	my $sequence = <IN> ;
	chomp $sequence;
	
	my $sequence_length = length($sequence);
	
	my $count;

	if ( index($name, "-")==-1 ) {
		#Assume this file is not collapsed, just count each seqeunce as 1
		$count = 1 ;
	} else  {
		#Assume file is collapsed (that is - sequence-ID has two numbers with a separating dash)
		(undef, $count) = $name =~ /^\>(\d+)\-(\d+)$/ ;

		# If the match failed, treat this fasta as not collapsed;
		$count = 1 if not defined $count ;
	}
	
	$histogram{$sequence_length} += $count ;
}

#Textual Output
if (0) {
	print "Length\tCount\n";
	foreach my $length_key ( sort { $a <=> $b } keys %histogram ) {
		print $length_key,"\t", $histogram{$length_key},"\n";
	}
	exit 0;
}

## Build the data as required by GD::Graph::bars.
## Data list has two items (each item is itself a list)
##   1. a list of x-axis labels (these are the keys from the histogram)
##   2. a list of values
my @data = (
	[ sort { $a <=> $b } keys %histogram ],
	[ map { $histogram{$_} } sort { $a <=> $b } keys %histogram ] ) ;

my $graph = new GD::Graph::bars (1000,800);

$graph->set(
	x_label => 'Length',
	y_label => 'Amount',
	title => 'Sequences lengths Distribution (after clipping)',
	bar_spacing => 10,
	transparent => 0,
	t_margin => 10,
	y_tick_number => 20,
	y_long_ticks => 1,
	)  or die $graph->error;

$graph->plot(\@data) or die $graph->error;
print OUT $graph->gd->png;

close IN;
close OUT;
