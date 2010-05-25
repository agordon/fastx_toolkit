#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

sub scan_upto_new_read()
{
	while ( my $line = <> ) {
		chomp $line ;
		return $line if ( $line =~ /^@/ );
	}
	return undef;
}

my $line;
while( 1 ) {
	# skip until the next read-id is found, if needed.
	if ( ! defined $line ) {
		$line = scan_upto_new_read();
		last unless $line; #end of file?
	}
	# Now, "$line" is a valid read-id.

	# Get the next line, should be nucleotides
	my $nuc = <>;
	last unless $nuc;
	chomp $nuc;

	if ($nuc !~ /^[ACGTN]+$/) {
		# Invalid nucleotide string - report and skip to the next read
		print STDERR "Skipping read '$line' (bad nucleotide string) at line $.\n";
		if ($nuc =~ /^@/ ) {
			#This line it is a new read by itself
			$line = $nuc ;
		} else {
			$line = undef ; # force a re-read until the next read-id is found
		}
		next;
	}

	# Get the next line, should be read-ID with "+" prefix
	my $read_id_2 = <>;
	last unless $read_id_2;
	chomp $read_id_2;

	if ($read_id_2 !~ /^\+/) {
		# Invalid read-id string - report and skip to the next read
		print STDERR "Skipping read '$line' (bad read-id-2 string) at line $.\n";
		if ($read_id_2 =~ /^@/ ) {
			#This line it is a new read by itself
			$line = $read_id_2;
		} else {
			#this is a random invalid line
			$line = undef ; # force a re-read until the next read-id is found
		}
		next;
	}

	# TODO:
	#   validate 'read-id-2' - it should be identical to the first read (except for the prefix character)

	# Get the next line, should be the quality scores
	my $quality_scores = <>;
	last unless $quality_scores;
	chomp $quality_scores ;

	if ( ! (
		( length($quality_scores) == length($nuc) )
		||
		( $quality_scores =~ /^[\d\- ]+$/ ) ) ) {
		# Invalid quality-scores string - report and skip to the next read
		print STDERR "Skipping read '$line' (bad quality scores) at line $.\n";
		if ($quality_scores =~ /^@/ ) {
			#This line it is a new read by itself
			$line = $quality_scores;
		} else {
			#this is a random invalid line
			$line = undef ; # force a re-read until the next read-id is found
		}
		next;
	}

	# TODO
	#   validate content of quality-scores
	if ( $quality_scores =~ /^[\d\- ]+$/ ) {
		my @scores = split /\s+/, $quality_scores;

		if ( scalar(@scores) != length($nuc) ) {
			print STDERR "Skipping read '$line' (bad numeric quality scores - wrong number of elements) at line $.\n";
			$line = undef;
			next
		}
	}

	##
	## If we got here - the read is OK, print it
	##
	print $line, "\n";
	print $nuc, "\n";
	print $read_id_2, "\n";
	print $quality_scores, "\n";

	$line = undef;

}
