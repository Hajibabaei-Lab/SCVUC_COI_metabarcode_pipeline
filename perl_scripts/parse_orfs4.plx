#!/path/to/miniconda3/envs/SCVUCv4.3/bin/perl
# Feb. 3/2020 by Teresita M. Porter
# Consolidate nt and aa ORFs for tranalign, no need to remove outliers
# Script to grab matching orfs from nt and aa files
# calculate cutoffs to remove nt orfs with length outliers
# output matching filtered nt and aa orfs
# USAGE perl parse_orfs.plx avg_orf_nt.fasta avg_orf_aa.fasta

use strict;
use warnings;
use Data::Dumper;

# declare var
my $otu;
my $orf;
my $length;
my $ntseq;
my $aaseq;
#my $lq = 25; #lower quartile
#my $uq = 75; #upper quartile
#my $percentile25;
#my $percentile75;
#my $iqr;
my $longest;
#my $min;
#my $max;
my $i=0;
my $outfile1 = $ARGV[0].".filtered";
my $outfile2 = $ARGV[1].".filtered";


# declare array
my @nt;
my @aa;
my @lengths;
my @otus;
my @longest;

# declare hash
my %ntLength; # key1 = otu, key2 = orf, value = length
my %ntSeq;
my %aaLength;
my %aaSeq;

my %match; # key1 = otu, key2 = orf, value = length

# read in files from ORFfinder
open (IN, "<", $ARGV[0]) || die "Can't open orf_nt.fasta: $!\n";
@nt = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Can't open orf_aa.fasta: $!\n";
@aa = <IN2>;
close IN2;

# hash length and seq from orf_nt.fasta
my ($ref1, $ref2) = parse_nt_orf(\@nt,\%ntLength, \%ntSeq);
%ntLength = %{$ref1};
%ntSeq = %{$ref2};

# hash length and seq from orf_aa.fasta
($ref1, $ref2) = parse_aa_orf(\@aa, \%aaLength, \%aaSeq);
%aaLength = %{$ref1};
%aaSeq = %{$ref2};

# keep records that match between nt and aa hashes
foreach $otu ( keys %ntLength) {
	foreach $orf (keys %{ $ntLength{$otu} }) {
		if (exists $aaLength{$otu}{$orf}) {
			$length = $ntLength{$otu}{$orf};
			$match{$otu}{$orf} = $length;
			push(@lengths, $length);
		}
	}
}

# calculate percentiles
#$percentile25 = get_percentile(\@lengths, $lq);
#$percentile75 = get_percentile(\@lengths, $uq);

# calculate interquartile range (IQR)
#$iqr = $percentile75 - $percentile25;

# calculate lower cutoff: 25th percentile - 1.5*IQR
#$min = $percentile25 - ($iqr * 1.5);

# calculate upper cutoff: 75th percentile + 1.5*IQR
#$max = $percentile75 + ($iqr * 1.5);

# create outfiles
open (OUT1, ">>", $outfile1) || die "Cannot open orf_nt.fasta.filtered:$!\n";
open (OUT2, ">>", $outfile2) || die "Cannot open orf_aa.fasta.filtered:$!\n"; 

foreach $otu ( keys %match ) {
#	foreach $orf ( keys %{ $match{$otu} } ) {
		# for each otu, sorth the orfs by ascending length
		@longest =  sort { $match{$otu}{$a} <=> $match{$otu}{$b} } keys %{ $match{$otu} };
		$longest = $longest[-1]; #inner key to longest is at bottom of array
		$length = $match{$otu}{$longest};

#		if ($length >= $min && $length <= $max) {
			$ntseq = $ntSeq{$otu}{$longest};	
			print OUT1 ">".$otu."\n".$ntseq."\n";
			$aaseq = $aaSeq{$otu}{$longest};
			print OUT2 ">".$otu."\n".$aaseq."\n";
#		}
#	}
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
my $orfaccession;
my @orfaccession;
my $orf;
my $accession;
my $orfbold;
#my $filter;
#my $PCRStep;
#my $amplicon;
#my $moleculeFilterPCRStepAmpliconOtu;

		@line = split(/ /, $line); #process headier line
			# 0 - >lcl|KR389058:1-588
			# 1 - ORF1_KR389058:0:588
		$orfline = $line[1];

		@orfline = split(/:/, $orfline);
			# 0 - ORF1_KR389058 or sometimes ORF1_BOLD:AAE3122:0:86
			# 1 - 0
			# 2 - 588

		if ($orfline =~ /BOLD/) {
			$orfbold = $orfline[0]; # ORF1_BOLD
			$accession = $orfline[1]; # AAE3122
			$start = $orfline[2]; # 0
			$stop = $orfline[3]; # 86
			$length = $stop - $start + 1;
			$orfaccession = $orfbold.":".$accession; # ORF1_BOLD:AAE3122
		}
		else {
			$orfaccession = $orfline[0];
			$start = $orfline[1];
			$stop = $orfline[2];
			$length = $stop - $start + 1;
		}

		@orfaccession = split(/_/, $orfaccession); # ORF1_KR389058 or ORF1_BOLD:AAE3122
			# 0 - ORF1
			# 1 - accession
		$orf = $orfaccession[0];
		$accession = $orfaccession[1];
#		$filter = $orfotu[2];
#		$PCRStep = $orfotu[3];
#		$amplicon = $orfotu[4];
#		$otu = $orfotu[5];
#		$moleculeFilterPCRStepAmpliconOtu = $molecule."_".$filter."_".$PCRStep."_".$amplicon."_".$otu;

		return ($orf, $accession, $length);
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
# create subroutine to parse FASTA header
sub parse_aa_header {

my $line = $_[0];
my @line;
my $orfline;
my @orfline;
my $start;
my $stop;
my $length;
my $orfaccession;
my @orfaccession;
my $orf;
my $accession;
my $orfbold;
#my $filter;
#my $PCRStep;
#my $amplicon;
#my $otu;
#my $moleculeFilterPCRStepAmpliconOtu;

		@line = split(/ /, $line); #process header line
		# 0 - >lcl|ORF1_KR389058:0:587 or sometimes >lcl|ORF1_BOLD:AAE3122:0:86
		# 1 - unnamed
		# 2 - protein
		# 3 - product,
		#4 - partial
		$orfline = $line[0];

		@orfline = split(/:/, $orfline);

		if ($orfline =~ /BOLD/) {
			# 0 - >lcl|ORF1_BOLD
			# 1 - AAE3122
			# 2 - 0
			# 3 - 86
			$orfbold = $orfline[0]; # lcl|ORF1_BOLD
			$orfbold =~ s/^>lcl\|//g; # ORF1_BOLD
			$accession = $orfline[1]; # AAE3122
			$orfaccession = $orfbold.":".$accession; # ORF1_BOLD:AAE3122
			$start = $orfline[2]; # 0
			$stop = $orfline[3]; # 86
			$length = $stop - $start + 1;
		}
		else {
			# 0 - >lcl|ORF1_KR389058
			# 1 - 0
			# 2 - 587
			$orfaccession = $orfline[0];
			$orfaccession =~ s/^>lcl\|//g;
			$start = $orfline[1];
			$stop = $orfline[2];
			$length = $stop - $start + 1;
		}

		@orfaccession = split(/_/, $orfaccession);
			# 0 - ORF1
			# 1 - accession
		$orf = $orfaccession[0];
		$accession = $orfaccession[1];
#		$filter = $orfotu[2];
#		$PCRStep = $orfotu[3];
#		$amplicon = $orfotu[4];
#		$otu = $orfotu[5];
#		$moleculeFilterPCRStepAmpliconOtu = $molecule."_".$filter."_".$PCRStep."_".$amplicon."_".$otu;

		return ($orf, $accession, $length);
 }

#######################################################
# parse out lengths and seqs from orf_aa.fasta

sub parse_aa_orf {

my @aa = @{$_[0]}; # access array reference
my %aaLength = %{$_[1]};
my %aaSeq = %{$_[2]};
my $i = 0;
my $line;
my $flag = 0;
my $orf;
my $otu;
my $length;
my $aaseq;

	while ($aa[$i]) {
		$line = $aa[$i];
		chomp $line;

		if ($flag==0 && $line =~/^>/) {
			($orf, $otu, $length) = parse_aa_header($line);
			$aaLength{$otu}{$orf} = $length;
			$flag = 1;
			$i++;
			next;
		}
		elsif ($flag==1) {
			$aaseq = $line;
			$flag = 2;
			$i++;
			next;
		}
		elsif ($flag==2 && $line !~ /^>/) {
			$aaseq = $aaseq.$line;
			$i++;
			next;
		}
		else {
			$flag = 0;
			$aaSeq{$otu}{$orf} = $aaseq;
			next;
		}
	}
	$i=0;

	# add last seq to hash
	$aaSeq{$otu}{$orf} = $aaseq;

	return(\%aaLength, \%aaSeq);

}

#######################################################
# Calculate percentile

#sub get_percentile {

#my @lengths = @{$_[0]};
#my $percentile = $_[1];

#my @sorted_lengths;
#my $array_length;
#my $decimal = $percentile/100;
#my $index;
#my $value;

#	@sorted_lengths = sort { $a <=> $b } @lengths;
#	$array_length = scalar(@sorted_lengths);
#	$index = int($array_length * $decimal);
#	$value = $sorted_lengths[$index];
#	return($value);

#}
