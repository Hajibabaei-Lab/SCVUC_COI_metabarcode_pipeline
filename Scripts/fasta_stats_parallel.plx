#!/usr/bin/perl
#June 28, 2018 by Teresita M. Porter
#Script calculates read stats from a FASTA-formatted file
#USAGE $perl fasta_stats_parallel.plx infile 

use strict;
use warnings;

use Statistics::Lite qw(:all);

#declare variables
my $filename;
my $base;
my $line;
my $flag=0;
my $seq1;
my $seq2;
my $seq;
my $i=0;
my $length;
my $min;
my $max;
my $mean;
my $mode;
my $j=0;

#declare array
my @seq;
my @split;
my @length;
my @filename;
my @fasta;

if (-z $ARGV[0] ) { #check whether filesize is 0
	print "File is empty\n";
}
else {

open (FASTA, "<", $ARGV[0]) || die "Error cannot open file: $!\n";
@fasta = <FASTA>;
close FASTA;

#customize filename parsing here if necessary
$filename = $ARGV[0];
chomp $filename;
#$base = substr $filename, 0, 8;
@filename = split(/\./,$filename);
$base = $filename[0];

while($fasta[$j]){
	$line = $fasta[$j];
	chomp $line;

	if ($flag==0){

		if ($line =~ /^>/){
			$j++;
			next;
		}
		elsif ($line =~ /^s+/) { ### skip blank lines ###
			$j++;
			next;
		}
		else {
			$seq1 = $line;
			$flag=1;
			$j++;
		}
		
	}	
	elsif ($flag==1){
		
		if ($line =~ /^>/){
			$flag=0;
			push (@seq, $seq1);
			$j++;
		}
		elsif ($line =~ /^s+/) { ### skip blank lines ###
			$j++;
			next;
		}
		else {
			$seq2 = $line;
			$seq1 = $seq1.$seq2;
			$j++;
		}
	}
}
$j=0;
push (@seq, $seq1);#don't forget to add last seq in file!
	
while ($seq[$i]){
	@split = split(//, $seq[$i]);
	$length = scalar(@split);
	push (@length, $length);
	@split =();#empty array
	$i++;
}

foreach (@length){
	my $element = $_;
}

$min = min (@length);
$max = max (@length);
$mean = mean (@length);
$mode = mode (@length);
my $num = scalar(@seq);

#Expected headers "Filename\tNumSeqs\tMin\tMax\tMean\tMode\n";
print $base."\t".$num."\t".$min."\t".$max."\t".$mean."\t".$mode."\n";
}
