#!/home/terri/miniconda3/envs/myenv/bin/perl
# Dec. 17, 2019 by Teresita M. Porter
# Accomodate FASTA header names ORF_otu only for the SCVUC pipeline
# calculate cutoffs to remove nt orfs with length outliers
# output filtered longest orfs (putative outliers removed)
# USAGE perl parse_orfs.plx cds.fasta.tmp cds.fasta

use strict;
use warnings;
use Data::Dumper;

# declare var
my $otu;
my $orf;
my $length;
my $ntseq;
my $lq = 25; #lower quartile
my $uq = 75; #upper quartile
my $percentile25;
my $percentile75;
my $iqr;
my $longest;
my $min;
my $max;
my $i=0;
my $outfile1 = $ARGV[1];

# declare array
my @nt;
my @lengths;
my @otus;
my @longest;

# declare hash
my %ntLength; # key1 = otu, key2 = orf, value = length
my %ntSeq;
my %ntLength_longest;
my %ntSeq_longest;

my %match; # key1 = otu, key2 = orf, value = length

# read in files from ORFfinder
open (IN, "<", $ARGV[0]) || die "Can't open orf_nt.fasta: $!\n";
@nt = <IN>;
close IN;

# hash length and seq from orf_nt.fasta
my ($ref1, $ref2) = parse_nt_orf(\@nt,\%ntLength, \%ntSeq);
%ntLength = %{$ref1};
%ntSeq = %{$ref2};

# get longest orfs
foreach $otu ( keys %ntLength ) {
	# ascending sort
	@longest =  sort { $ntLength{$otu}{$a} <=> $ntLength{$otu}{$b} } keys %{ $ntLength{$otu} };
	$longest = $longest[-1]; #inner key to longest is at bottom of array
	# get length and seq from original hashes
	$length = $ntLength{$otu}{$longest};
	$ntseq = $ntSeq{$otu}{$longest};
	# put longest orf length and seq into new hashes 
	$ntLength_longest{$otu}{$longest} = $length;
	$ntSeq_longest{$otu}{$longest} = $ntseq;
}

# grab length of all longest orfs
foreach $otu ( keys %ntLength_longest) {
	foreach $orf (keys %{ $ntLength_longest{$otu} }) {
		$length = $ntLength_longest{$otu}{$orf};
		push(@lengths, $length);
	}
}

# calculate percentiles
$percentile25 = get_percentile(\@lengths, $lq);
$percentile75 = get_percentile(\@lengths, $uq);

# calculate interquartile range (IQR)
$iqr = $percentile75 - $percentile25;

# calculate lower cutoff: 25th percentile - 1.5*IQR
$min = $percentile25 - ($iqr * 1.5);

# calculate upper cutoff: 75th percentile + 1.5*IQR
$max = $percentile75 + ($iqr * 1.5);

# testing
#print "min $min\t max $max\n";

# create outfiles
open (OUT1, ">>", $outfile1) || die "Cannot open orf_nt.fasta.filtered:$!\n";

foreach $otu ( keys %ntLength_longest ) {
	foreach $orf (keys %{ $ntLength_longest{$otu} }) {
		$length = $ntLength_longest{$otu}{$orf};

		if ($length >= $min && $length <= $max) {
			$ntseq = $ntSeq_longest{$otu}{$orf};	
			print OUT1 ">".$otu."\n".$ntseq."\n";
		}
	}
}

#######################################################
# create subroutine to parse FASTA header
sub parse_nt_header {

my $line = $_[0];
my @line;
my $orfline;
my @orfline;
my $start;
my $stop;
my $orfotu;
my @orfotu;
my $orf;
my $otu;

		@line = split(/ /, $line); #process headier line
			# 0 - >lcl|Otu1:2:310
			# 1 - ORF1_Otu1:1:309
		$orfline = $line[1];

		@orfline = split(/:/, $orfline);
			# 0 - ORF1_Otu1
			# 1 - 1
			# 2 - 309
		$start = $orfline[1];
		$stop = $orfline[2];
		$length = $stop - $start + 1;
		$orfotu = $orfline[0];

		@orfotu = split(/_/, $orfotu);
			# 0 - ORF1
			# 1 - Otu1
		$orf = $orfotu[0];
		$otu = $orfotu[1];

		return ($orf, $otu, $length);
 }

#######################################################
# parse out lengths and seqs from orf_nt.fasta

sub parse_nt_orf {

my @nt = @{$_[0]}; # access array reference
my %ntLength = %{$_[1]};
my %ntSeq = %{$_[2]};
my $i = 0;
my $line;
my $flag = 0;
my $orf;
my $otu;
my $length;
my $ntseq;

	while ($nt[$i]) {
		$line = $nt[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			($orf, $otu, $length) = parse_nt_header($line);
			$ntLength{$otu}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1) {
			$ntseq = $line;
			$flag = 2;
			$i++;
			next;
		}
		elsif ($flag==2 && $line !~ /^>/) {
			$ntseq = $ntseq.$line;
			$i++;
			next;
		}
		else {
			$flag = 0;
			$ntSeq{$otu}{$orf} = $ntseq;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$ntSeq{$otu}{$orf} = $ntseq;

	return(\%ntLength, \%ntSeq);

}

#######################################################
# Calculate percentile

sub get_percentile {

my @lengths = @{$_[0]};
my $percentile = $_[1];

my @sorted_lengths;
my $array_length;
my $decimal = $percentile/100;
my $index;
my $value;

	@sorted_lengths = sort { $a <=> $b } @lengths;
	$array_length = scalar(@sorted_lengths);
	$index = int($array_length * $decimal);
	$value = $sorted_lengths[$index];
	return($value);

}
