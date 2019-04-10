#!/usr/bin/perl

use warnings;
use strict;

my %genehash;

if (defined($ARGV[0])) {
	print STDERR "Usage: remove_graph_duplicates.pl <GRAPH\n\nChecks a given graph input for duplicates and removes them (species only)\nReads from STDIN\n";
	exit;
}

my $a = "unk";
my $b = "unk";
while (<STDIN>) {
	my @row = split(/\s+/);
	if ($_ =~ /^#/) {
		unless ($row[1]) {die();}
		unless ($row[1] eq "file_a" || $row[1] eq "a") {
			$a = $row[1];
			$b = $row[2];
		}
	}
	else {
		unless ($row[1]) {next;}

#		print "$a: $row[0]\n$b: $row[1]\n";

		if (defined($genehash{$row[0]})) {
			if ($genehash{$row[0]} ne $a) {
			print "$row[0] was used with $genehash{$row[0]} and $a\n";
			}
		}
		else {
			$genehash{$row[0]} = $a;
		}
		if (defined($genehash{$row[1]})) {
			if ($genehash{$row[1]} ne $b) {
			print "$row[1] was used with $genehash{$row[1]} and $b\n";
			}
		}
		else {
			$genehash{$row[1]} = $b;
		}
	}
}
