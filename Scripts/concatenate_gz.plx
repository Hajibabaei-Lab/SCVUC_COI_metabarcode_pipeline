#!/usr/bin/perl
#Nov. 23/17 by Teresita M. Porter
#Script to combine sequence files when a run has been repeated
#USAGE perl concatenate_gz.plx

use strict;
use warnings;

#var
my $dir;
my $i=0;
my $file;
my $sample;
my $list;
my $newlist;
my $value;
my $filename;
my $R;
my $newfilename; ## see hardcoded GRDI formatting below L125
my $newfilename2; ## see hardcoded GRDI formatting below L137
my $dirfilename;

#array
my @files;
my @index;
my @list;
my @filename;

#hash
my %filename; #key=sample, value=filename
my %sample; #key=sample, value=1
my %R1; #key=sample, value=filename
my %R2;

print "Enter name of directory containing sequences (including final '/':\n";

$dir = <STDIN>;
chomp $dir;

opendir (DIR, $dir) || die "Can't open directory: $!\n";
@files = readdir DIR;

#hash of sample counts indexed by samples
#hash of filenames indexed by samples
while ($files[$i]) {
	$file = $files[$i];
	chomp $file;
	#print $file."\n";
	
	if ($file =~ /^\./) {
		$i++;
		next;
	}
	
	if ($file !~ /LV/) {
		$i++;
		next;
	}

	if ($file =~ /ECOSM/) {
		@index = split("_",$file);
		$sample = $index[0];
	}
	else {
		@index = split("_",$file);
		$sample = $index[1];
	}

	if (exists $filename{$sample}) {
		$list = $filename{$sample};
		$newlist = $list."|".$file;
		$filename{$sample} = $newlist;
	}
	else {
		$filename{$sample} = $file;
	}
	
	if (exists $sample{$sample}) {
		$value = $sample{$sample};
		$value++;
		$sample{$sample} = $value;
	}
	else {
		$sample{$sample} = 1;
	}	
	
	$i++;
}
$i=0;

while (($sample, $value) = each(%sample)) {
	if ($value > 2) {
		$list = $filename{$sample};
		@list = split(/\|/,$list);
		
		foreach $filename (@list){
			@filename = split('_', $filename);
			$R = $filename[4];

			if ($R =~ /R1/) {
				if (exists $R1{$sample}) {
					$list = $R1{$sample};
					$newlist = $list."|".$filename;
					$R1{$sample} = $newlist;
				}
				else {
					$R1{$sample} = $filename;
				}
			}
			elsif ($R =~ /R2/) {
				if (exists $R2{$sample}) {
					$list = $R2{$sample};
					$newlist = $list."|".$filename;
					$R2{$sample} = $newlist;
				}
				else {
					$R2{$sample} = $filename;
				}
			}
			else {
				print "Error with R: $R\n";
			}
		}
		
		$newfilename = $dir."GRDI-ECO_".$sample."_F230_COMB_R1.fq.gz";
		open (CAT, ">>", $newfilename) || die "Error cannot open $newfilename: $!\n";
		close CAT;

		$list = $R1{$sample};
		@list = split(/\|/,$list);
		foreach $filename (@list) {
			$dirfilename = $dir.$filename;
			system("cat $dirfilename >> $newfilename");
			unlink $dirfilename;		
		}

		$newfilename2 = $dir."GRDI-ECO_".$sample."_F230_COMB_R2.fq.gz";
		open (CAT, ">>", $newfilename2) || die "Error cannot open $newfilename2: $!\n";
		close CAT;

		$list = $R2{$sample};
		@list = split(/\|/,$list);

		foreach $filename (@list) {
			$dirfilename = $dir.$filename;
			system ("cat $dirfilename >> $newfilename2");
			unlink $dirfilename;
		}
	}
}
