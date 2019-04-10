#!/usr/bin/perl

use warnings;
use strict;

print "graph G {\n";
while (<STDIN>) {
	chomp;
	if ($_ =~ /#/) {next;}
	my @x = split(/\s/);
	if ($x[1]) {
		print " $x[0] -- $x[1];\n"
	}
}
print " }";
