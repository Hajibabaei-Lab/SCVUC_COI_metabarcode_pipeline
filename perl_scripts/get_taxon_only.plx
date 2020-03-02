#!/path/to/miniconda3/envs/SCVUCv4.3/bin/perl
# March 2, 2020 by Teresita M. Porter
# Script to filter cat.denoised.nonchimeras to previously chosen taxon only (see config file)
# USAGE perl get_taxon_only.plx taxon.zotus cat.denoised.nonchimeras > cat.denoised.nonchimeras.taxon

use strict;
use warnings;
use Data::Dumper;

# vars
my $i=0;
my $line;
my $zotu;
my $oldseq="";
my $seq;
my $flag=0;

# arrays
my @tax;
my @fas;

# hashes
my %fasta; #key = zotu, val = seq

open (TAXON, "<", $ARGV[0]) || die "Cannot open list of arthropod zotus: $!\n";
@tax = <TAXON>;
close TAXON;

open (FAS, "<", $ARGV[1]) || die "Cannot open cat.denoised.nonchimeras FASTA file: $!\n";
@fas = <FAS>;
close FAS;

# hash fasta file
while ($fas[$i]){
	$line = $fas[$i];
	chomp $line;

	if ($line =~ /^>/) {
		if ($flag == 0) {
			$line =~ s/^>//;
			$zotu = $line;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag == 1) {
			$flag = 0;
			$fasta{$zotu} = $oldseq;
			$oldseq="";
			next;
		}
	}
	elsif ($line !~ /^>/ && $flag == 1) {
		$seq = $oldseq.$line;
		$oldseq = $seq;
		$i++;
		next;
	}
}
$i=0;

# add final record
$fasta{$zotu} = $oldseq;

while ($tax[$i]) {
	$zotu = $tax[$i];
	chomp $zotu;

	if (exists $fasta{$zotu}) {
		$seq = $fasta{$zotu};
		print ">$zotu\n$seq\n";
	}
	$i++;
}
$i=0;
