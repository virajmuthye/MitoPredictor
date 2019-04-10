#!/usr/bin/perl

use warnings;
use strict;

my %filehash;
my $flag;

if (defined($ARGV[0])) {
	print STDERR "Usage: remove_graph_duplicates.pl <GRAPH\n\nChecks a given graph input for duplicates and removes them (species only)\nReads from STDIN\n";
	exit;
}

while (<STDIN>) {
	if ($_ =~ /^#/) {
		$flag = 1;
		my @row = split(/\s+/);
		unless ($row[1]) {die();}
		unless ($row[1] eq "file_a" || $row[1] eq "a") {
			my $a = $row[1];
			my $b = $row[2];
			my $id = join(",",sort($a,$b));

			if (++$filehash{$id} > 1) {
				$flag = 0;
				print STDERR "Found duplicate: $a vs $b\n";
			}
		}
	}
	if ($flag) {
		print $_;
	}
}
