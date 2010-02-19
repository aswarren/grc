#!/usr/bin/perl -w
use Getopt::Std;
#the following is a perl script to parse a fasta file based on user choice/fasta headers
#$YorN=1; #boolean for yes or no answer


open OPfile, "> ".$ARGV[-1];

$TotalCount=0;
print "Merging fasta files\n";

foreach $a (@ARGV){#open for loop
	if($a eq $ARGV[-1]){
		last;
	}
	$InHandle=$a;#set the input file name
	open $IPfile, "< $InHandle";
	
	my $count=0;
	while (<$IPfile>) {
		chomp($_);
		print OPfile "$_\n"; #print line
		if (/>/){ #header
			$count++;}
	} #close while loop
	
	
	print "$InHandle: $count records\n";
	$TotalCount=$count+$TotalCount;
	
	close $IPfile;

}#close for loop

	
	
	
	
	close OPfile;
	
	print "TotalCount: $TotalCount \n\n";



