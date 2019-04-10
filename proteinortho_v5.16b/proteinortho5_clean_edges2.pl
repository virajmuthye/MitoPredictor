#!/usr/bin/perl

use strict;
use warnings "all";

if (!defined($ARGV[1])) {
	print STDERR "proteinortho5_clean_edges.pl RM_LIST EDGELIST(S) >CLUSTERED_EDGELIST \nReads Proteinortho outfile and its initial edge list and calculates\nthe remaining edge list after the step of clustering.\n";
	exit;
}


# Proteinortho Output
print STDERR "Cleaning edge list...\n";
# Reading removal list
my %map;
open(IN,"<$ARGV[0]") || die("Could not open file $ARGV[0]: $!");
while (<IN>) {
	chomp;
	my ($ta, $sa, $tb, $sb) = split("\t",$_,4);	# name, species, name, species
	unless ($sb) {die("Line does not match filter list pattern\n");}

	$ta .= " ".$sa;
	$tb .= " ".$sb;
	my ($a, $b) = sort($ta,$tb); # bidirectional edges store once is fine
	$map{$a.' '.$b} = 1;
	
}
close(IN);

# Edgelist
shift(@ARGV);
my $rm = 0;
my $all = 0;
# For all blastgraphs
my $filea = "";
my $fileb = "";
foreach my $file (@ARGV) {
	open(IN,"<$file") || die("Could not open file $file: $!");
	while (<IN>) {
		if ($_ =~ /^#/) {
			print $_; 
			chomp;
			$_ =~ s/^#\s*//;
			my ($a,$b,undef) = split("\t",$_,3);
			if ($a eq "a" || $a eq "file_a") {next;}
			unless (defined($b)) {die("Line does not match Proteinortho graph format. Filename missing.\n");}
			# Keep track of files/species
			$filea = $a;
			$fileb = $b;
		}
		else {
			unless ($_ =~ /\t/) {next;}
			my ($ta,$tb,undef) = split("\t",$_,3);
			unless ($tb) {die("Line does not match Proteinortho graph format\n");}
			$all++;
	
			$ta .= " ".$filea;
			$tb .= " ".$fileb;

			my ($a, $b) = sort($ta,$tb); # bidirectional edges store once is fine
			($a,$b) = sort($a, $b);
			if (exists($map{$a.' '.$b})) {
				$rm++;
				next;
			}
	
	
#			if ($a eq $b) {$all--; next;}	# Selfhits are ignored
			print $_;
		}
	}
	close(IN);
}

print STDERR "Removed $rm / $all edges\n";
print STDERR "Done.\n";

## TODO: Adjust to remove selfblast hits and watch
