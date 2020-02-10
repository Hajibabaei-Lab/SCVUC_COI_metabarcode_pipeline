#!/path/to/miniconda3/envs/myenv.3/bin/perl
# Feb. 10, 2020 by Teresita M. Porter
# Script to add longest arthropod ORF sequences to the result file
# expect a STRICT fasta file
# USAGE perl add_longest_orfs_to_csv.plx longest.orfs.fasta rdp.csv.tmp > rdp.csv.tmp2

use strict;
use warnings;
use Data::Dumper;

# vars
my $i=0;
my $line;
my $zotu;
my $j;
my $seq;
my $add;

# arrays
my @orfs;
my @tax;
my @line;

# hashes
my %orfs; #key = accession, value = seq

open (ORFS, "<", $ARGV[0]) || die "Cannot open longest orfs fasta file:$!\n";
@orfs = <ORFS>;
close ORFS;

open (TAX, "<", $ARGV[1]) || die "Cannot open rdp taxonomy file: $!\n";
@tax = <TAX>;
close TAX;

# hash orfs
while ($orfs[$i]) {
	$line = $orfs[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$zotu = $line;
		$j = $i+1;
		$seq = $orfs[$j];
		chomp $seq;
		$orfs{$zotu} = $seq;
	}
	$i+=2;
	
}
$i=0;

#print Dumper(\%orfs); # tesing

# hash tax
while ($tax[$i]) {
	$line = $tax[$i];
	chomp $line;

	if ($i==0) {
	print "GlobalESV,LongestORFseq,SampleName,ESVsize,Strand,Root,RootRank,rBP,SuperKingdom,SuperKingdomRank,skBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP\n";	
		$i++;
		next;
	}

	@line = split(/,/, $line);
	$zotu = shift(@line);

	if (exists $orfs{$zotu} ) {
		$seq = $orfs{$zotu};
		$add = $zotu.",".$seq;
		unshift @line, $add;
		$line = join ',', @line;
		print $line."\n"; #tesing
	}
	$i++;
}
$i=0;
