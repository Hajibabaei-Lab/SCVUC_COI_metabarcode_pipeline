#!/path/to/miniconda3/envs/myenv/bin/perl
# Nov. 7, 2019 by Teresita M. Porter
# Script to grab longest CDS from ORFfinder nucleotide output
# Keep CDS with lengths within 25th percentile + 1.5*IQR to 75th percentile + 1.5*IQR
# USAGE perl get_longest_orf.plx cds.fasta.tmp limits.txt cds.pdf cds.fasta.tmp2 cds.fasta

use strict;
use warnings;

# declare var
my $i=0;
my $line;
my $orf;
my $otu;
my $length;
my $length2;
my $flag=0;
my $seq;
my $seq1;
my $seq2;
my $longestOrf;
my $longestOrf2;
my $longest;
my $lower;
my $upper;
my $num;
my $num2;

my $limits = $ARGV[1]; # matches range in snakefile
my $plot = $ARGV[2]; # matches plot in snakefile
my $cds_out = $ARGV[3]; # matches cds_out in snakefile
my $cds_out2 = $ARGV[4]; # matches cds_out2 in snakefile

# declare array
my @in;
my @line;
my @otus;
my @orfs;
my @longestOrfs;
my @longestOrfs2;
my @orfs2;
my @sorted;
my @sorted2;
my @range;

# declare hash
my %length;
my %seq;

open (IN, "<", $ARGV[0]) || die "Can't open infile: $!\n";
@in=<IN>;
close IN;

# parse file to get fasta header stats first
while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($flag==0 && $line =~/^>/) {
		($orf, $otu, $length) = parse_header($line);
		$length{$otu}{$orf} = $length; #create hash of hashes
		$flag = 1;
		$i++;
		next;
	}
	elsif ($flag==1) {
		$seq = $line;
		$flag = 2;
		$i++;
		next;
	}
	elsif ($flag==2 && $line !~ /^>/) {
		$seq = $seq.$line;
		$i++;
		next;
	}
	else {
		$flag = 1;
		$seq{$otu}{$orf} = $seq;
		($orf, $otu, $length) = parse_header($line);
		$length{$otu}{$orf} = $length;
		$i++;
		next;
	}
}
$i=0;
		
# add last seq to hash
$seq{$otu}{$orf} = $seq;

# loop through each otu, get orf id for longest length, print otu and seq to outfile

# get unique otu keys
@otus = keys %length;

# open outfile for longest nt ORFs (not filtered, contains outliers)
open (OUT, ">>", $cds_out) || die "Error cannot open cds_out: $!\n";

# print longest ORFs before and after removing outliers
while ($otus[$i]) {
	$otu = $otus[$i];

	# get length of longest orfs
	@orfs =  sort { $length{$otu}{$a} <=> $length{$otu}{$b} } keys %{$length{$otu}};
	$longestOrf = $orfs[-1]; #inner key to longest is at bottom of array
	$length = $length{$otu}{$longestOrf};

	# get array with all longest orfs (there may be a few equally long orfs)
	foreach $orf (keys %{$length{$otu}}) {
		if ($length{$otu}{$orf} == $length) {
			push(@longestOrfs, $orf);
		}
		else {
			next;
		}
	}

	# get number of equally long orfs 
	$num = scalar(@longestOrfs);

	# ascending sort longest OrfsIDs to get the first one with longest length
	if ($num > 1) {
		@sorted = sort { number_strip($a) <=> number_strip($b) } @longestOrfs;
		$longestOrf = $sorted[0];
	}

	$seq1 = $seq{$otu}{$longestOrf};
	print OUT ">$otu\n$seq1\n";

	$i++;
	@longestOrfs=();
	next;
}
$i=0;

close OUT;

# run R script to plot hisogram and get cutoff lengths
# paths relative to original snakefile
system("Rscript R_scripts/length_hist.R $cds_out $plot $limits");
wait;

# get 1.5*IQR
open (IN2, "<", $limits) || die "Can't open limits.txt:$!\n";
@range = <IN2>;
close IN2;

$lower = $range[0];
chomp $lower;
$upper = $range[1];
chomp $upper;

# open outfile for longest nt ORFs, filtered, outliers removed
open(OUT2, ">>", $cds_out2) || die "Error cannot open cds_out2: $!\n";

# print longest ORFs before and after removing outliers
while ($otus[$i]) {
	$otu = $otus[$i];

	# get length of longest orfs
	@orfs2 =  sort { $length{$otu}{$a} <=> $length{$otu}{$b} } keys %{$length{$otu}};
	$longestOrf2 = $orfs2[-1]; #inner key to longest is at bottom of array
	$length2 = $length{$otu}{$longestOrf2};

	# get array with all longest orfs (there may be a few equally long orfs)
	foreach $orf (keys %{$length{$otu}}) {
		if ($length{$otu}{$orf} == $length2) {
			push(@longestOrfs2, $orf);
		}
	}

	$num2 = scalar(@longestOrfs2);

	# ascending sort longest OrfsIDs to get the first one with longest length
	if ($num2 > 1 ) {
		@sorted2 = sort { number_strip($a) <=> number_strip($b) } @longestOrfs2;
		$longestOrf2 = $sorted2[0];
	}

	if ($length2 >= $lower && $length2 <= $upper) {
		$seq2 = $seq{$otu}{$longestOrf2};
		print OUT2 ">$otu\n$seq2\n";
		$i++;
		next;
	}
	else {
		$i++;
		next;
	}
}
$i=0;

close OUT2;

#######################################################
# subroutine to get number from OrfID

sub number_strip {
    my $line = shift;
    my ($num) = $line =~ /(\d+)/;
    return $num;
}

#######################################################
# subroutine to parse FASTA header
sub parse_header {

my $line = $_[0];
my @line;
my $orfline;
my @orfline;
my $start;
my $stop;
my $length;
my $orfotu;
my @orfotu;
my $orf;
my $otu;

		@line = split(/ /, $line); #process header line
			# 0 - >lcl|Otu1:2-310
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

