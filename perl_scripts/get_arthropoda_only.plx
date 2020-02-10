#!/path/to/miniconda3/envs/myenv.3/bin/perl
# Feb. 10, 2020 by Teresita M. Porter
# Script to filter cat.denoised.nonchimeras to arthropoda only
# USAGE perl get_arthropoda_only.plx arthropoda.zotus cat.denoised.nonchimeras > cat.denoised.nonchimeras.arthropoda

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
my @arth;
my @fas;

# hashes
my %fasta; #key = zotu, val = seq

open (ARTH, "<", $ARGV[0]) || die "Cannot open list of arthropod zotus: $!\n";
@arth = <ARTH>;
close ARTH;

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

#print Dumper(\%fasta); # testing

while ($arth[$i]) {
	$zotu = $arth[$i];
	chomp $zotu;

	if (exists $fasta{$zotu}) {
		$seq = $fasta{$zotu};
		print ">$zotu\n$seq\n";
	}
	$i++;
}
$i=0;
