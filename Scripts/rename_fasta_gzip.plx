#!/usr/bin/perl
#June 28, 2018 by Teresita M. Porter
#Script to add sampe name from files with GRDI-formatting to the fasta header prior to global concatenation and denoising otherwise OTU tables aren't generated properly
#USAGE perl rename_fasta_gzip.plx file.fasta.gz

use strict;
use warnings;

#declare vary
my $filename;
my $base;
my $outfile;
my $line;
my $i=0;
my $newline;
my $temp;
my $flag=0;

#declare array
my @in;
my @filename;

if ($ARGV[0] =~ /gz$/) {
	open (IN, "gunzip -c $ARGV[0] |") || die "Cannot open gzipped infile: $!\n";
	@in=<IN>;
	close IN;
	$flag=1;
}
else {
	open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
	@in=<IN>;
	close IN;
}

$filename = $ARGV[0];
@filename = split(/\./,$filename); ##### edit delimiter type here #####
$base = $filename[0];

$outfile = $base.".renamed.fasta.gz";

open (OUT, "| gzip -c >> $outfile") || die "Cannot open gzipped outfile: $!\n";

while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g; 
		$newline = $base.";".$line;
		print OUT ">".$newline."\n";
	}
	else {
		print OUT "$line\n";
	}
	$newline=();
	$i++;
}
$i=0;

close OUT;
