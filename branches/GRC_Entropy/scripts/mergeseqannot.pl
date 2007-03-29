#!/usr/bin/perl -w
use Getopt::Std;
#the following is a perl script to merge AA sequence and annotation files 
#These must be provided in (annotation, sequence) pairs e.g. "./mergeseqannot.pl a1.faa a1.ptt a2.faa a2.ptt" with the last parameter being the output file
#e.g. "./mergeseqannot.pl a1.fasta a1.goa a2.fasta a2.goa"



my %Annotation =(); #hash to store annotation information
my $NumParam=@ARGV;#get the number of parameters
$NumParam=$NumParam-1;#subtract off output parameter
$TotalCount=0;
print "Merging sequence and annotation files\n";

open ($OPfile, "> ".$ARGV[-1])
	or die "Couldn't open output file for writing: $!\n";

if($ARGV[1]=~/.ptt/){
	$a=0;
	while($a+1<$NumParam && $ARGV[$a+1]=~/ptt$/i){#open parameter loop
		$InHandle=$ARGV[$a+1];#set the input file name
		open ($IPfile, "< $InHandle")
			or die "Couldn't open input file: $!\n";#open ptt file
		local @Lines=<$IPfile>;#get contents
		chomp(@Lines);
	
		
		$count=3;
		#Create hash
		while ($count<@Lines) {#for each line of the ptt (starting at 4th line)
			local @Terms=split(/\t/, $Lines[$count]);
			$Annotation{$Terms[3]}=$Lines[$count]; #hash the ptt line to the gi number
			$count++;
		} #close while loop
		
		
		print "$InHandle: $count records\n";
		$TotalCount=$count+$TotalCount;
		close $IPfile;
		@Lines=();

		if($a<$NumParam && $ARGV[$a]=~/faa$/i){
			$InHandle=$ARGV[$a];#set the input file name
			open ($IPfile, "< $InHandle")
				or die "Couldn't open input file: $!\n";#open faa file
			while(<$IPfile>){#open faa for loop
				if($_=~/>/){#if aa record header
					local @Terms=split(/\|+/, $_);#break up line to get ID
					local $Info=$Annotation{$Terms[1]};#use 2nd element (should be ID) to get ptt line
					if(defined($Info)){#if its found
						local @AnnTerms=split(/\t/,$Info);#split by tabs
						chomp @AnnTerms;
						$Function="$AnnTerms[-1] $AnnTerms[4]";
					}#close if in hash
					else{#else not in hash
						$Function=$Terms[-1];#everything after last | in faa header
					}
					$Function=~s/\|+|\-+/ /g; #replace | or - with a space
					$Function=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
					$Function=~ tr/A-Z/a-z/; #convert everything to lower case
					$Function=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
					$Function=~ s/^\s+//gm; #remove leading whitespace
					$Function=~ s/\s+$//; #remove trailing whitespace
					print $OPfile ">".$Terms[1]."\t$Function\t$InHandle\n";#print ">ID ProductDescription FileName"
				}#close if header
				else{#if its not a header
					chomp;
					print $OPfile $_."\n";
				}
			}#close while faa open
		}

		else{#else annotation sequence file mismatch
			die "error combining files\n";
		}
		
		%Annotation=();#clear hash
		$a=$a+2;#next pair
	}#close parameter loop
	
	print "TotalCount: $TotalCount \n\n";
}#close if ptt



elsif($ARGV[1]=~/.goa/){
	$IDPosition=1;
	$GOPosition=4;
	$ECodePosition=6;
	$a=0;
	while($a+1<$NumParam && $ARGV[$a+1]=~/goa$/i){#open parameter loop
		$InHandle=$ARGV[$a+1];#set the input file name
		open ($IPfile, "< $InHandle")
			or die "Couldn't open input file: $!\n";#open goa file
		local @Lines=<$IPfile>;#get contents
		chomp(@Lines);
		
		$count=0;
		#Create hash
		while ($count<@Lines) {#for each line of the goa file
			local @Terms=split(/\t/, $Lines[$count]);
			#print "$Terms[$IDPosition] $Terms[$GOPosition] $Terms[$ECodePosition] \n";
			if(defined $Annotation{$Terms[$IDPosition]}){#if this ID is in the hash
				$Annotation{$Terms[$IDPosition]}="$Annotation{$Terms[$IDPosition]} $Terms[$GOPosition] $Terms[$ECodePosition]"; #hash the goa line to the ID
			}
			else{
				$Annotation{$Terms[$IDPosition]}="";#clear the string
				$Annotation{$Terms[$IDPosition]}="$Terms[$GOPosition] $Terms[$ECodePosition]";
			}
			
			$count++;
		} #close while loop
		
		
		print "$InHandle: $count records\t";
		$TotalCount=$count+$TotalCount;
		close $IPfile;
		@Lines=();

		if($a<$NumParam && $ARGV[$a]=~/fasta$/i){
			$InHandle=$ARGV[$a];#set the input file name
			open ($IPfile, "< $InHandle")
				or die "Couldn't open input file: $!\n";#open fasta file
			while(<$IPfile>){#open faa for loop
				if($_=~/>/){#if aa record header
					local @Terms=split(/\s+/, $_);#break up line on whitespace to get ID
					local $ID=$Terms[0];
					$ID=~s/>//;#remove '>'
					if($ID=~/\|/){#if the ID contains a seperator
						@IDSplit=split(/\|/,$ID);
						$SecondID=$IDSplit[1];
						$ID=~s/\|$SecondID//g;#remove everything after seperator
						#print "$ID\n";
					}
					
					$OtherFunc=$_;
					$OtherFunc=~s/>$ID//; #remove ID from the rest of the annotation
					$OtherFunc=~s/\|+|\-+/ /g; #replace | or - with a space
					$OtherFunc=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
					$OtherFunc=~ tr/A-Z/a-z/; #convert everything to lower case
					$OtherFunc=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
					$OtherFunc=~ s/^\s+//gm; #remove leading whitespace
					$OtherFunc=~ s/\s+$//; #remove trailing whitespace
					local $Info=$Annotation{$ID};#use ID to get GO annotations
					if(defined($Info)){#if its found
						#local @AnnTerms=split(/\t+/,$Info);#split by tabs
						#chomp @AnnTerms;
						print $OPfile ">".$ID."\t$Info $OtherFunc\t$InHandle\n";#print ">ID ProductDescription FileName"
					}#close if in hash
					else{#else not in hash
						print $OPfile ">".$ID."\t$OtherFunc\t$InHandle\n";
					}
				}#close if header
				else{#if its not a header
					chomp;
					print $OPfile $_."\n";
				}
			}#close while fasta open
		}

		else{#else annotation sequence file mismatch
			die "error combining files\n";
		}
		
		%Annotation=();#clear hash
		$a=$a+2;#next pair
	}#close parameter loop
	
	print "TotalCount: $TotalCount \n\n";
}#close if goa

	
close OPfile;





