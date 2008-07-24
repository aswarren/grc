#!/usr/bin/perl -w
use Getopt::Std;
use Text::Balanced qw (extract_delimited);
#the following is a perl script to parse a ptt or goa/CP files for use as a reference by grc
#Assumes for example.goa there is example.CP in the same directory

getopt('i');# get the reference file


#open annotation
#open($Annotes, "< ./$opt_a")
#	or die "Couldn't open annotation for reading: $!\n";
my %Annotation =(); #hash to store annotation information

#open map
open($Map, "< $opt_i")
	or die "Couldn't open map for reading: $!\n";
#$Output="RefParse.txt";
#$Output=$opt_i;
#$Output=~s/.ptt/.parseptt/;

#open OPfile, "> $Output";
if($opt_i=~/.ptt/){
	local @Lines=<$Map>;#get contents
	chomp(@Lines);
	$count=3;
	
	while ($count<@Lines) {#for each line starting at the fourth line

		$TheLine=$Lines[$count];
		@Words=split(/\t/, $TheLine); #changed from \t+ to \t 10/08/06
		$Orientation=$Words[1];#the orientation effects the ordering of coordinates
		$ID=$Words[3];#fourth one in ptt is ID
		$ID=~ s/^\s+//gm; #remove leading whitespace
		$ID=~ s/\s+$//; #remove trailing whitespace
		#print $ID."\n";
		#$Coord=<$Map>;#get the genome description line
		#chomp($Coord=<$Map>);#get the coordinates line
		#@Coord=~ s/^\s+//gm; #remove leading whitespace
		$Function=$Words[-1];#annotation
		if($Words[4] ne "-"){#if the synonym exists
			$Function=$Function." $Words[4]";#add the synonym to the end
		}
		$Function=~s/\|+|\-+/ /g; #replace | or - with a space
		$Function=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
		$Function=~ tr/A-Z/a-z/; #convert everything to lower case
		$Function=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
		$Function=~ s/^\s+//gm; #remove leading whitespace
		$Function=~ s/\s+$//; #remove trailing whitespace
		
		@Cds=split(/\.\./, $Words[0]); #split the first element up into coordinates
		
		print "$ID\t";
		
		foreach $Cd (@Cds){
			$Cd=~ s/^\s+//gm; #remove leading whitespace
			$Cd=~ s/\s+$//; #remove trailing whitespace
			
		}#close foreach
		if($Orientation eq "+"){
			print "$Cds[0]\t";
			print "$Cds[1]\t";
		}
		else {
			print "$Cds[1]\t";
			print "$Cds[0]\t";
		}
			
		print $Function."\n";
		$count++;
	}#close while
	close $Map;
}#close ptt

#This section supports the merging of goa and .CP files for use as a reference
#Assumes for example.goa there is example.CP in the same directory
#if the Uniprot ID in goa annotation is not in the cooresponding CP file then the protein will have NotFound as the annotation
if($opt_i=~/.CP/){
	$IDColumngoa=1;
	$GOColumn=4;
	$ECodeColumn=6;
	$GoaFile=$opt_i;#create chromosome table name
	$GoaFile=~s/.CP/.goa/;#create chromosome table name

	local @Lines=<$Map>;#get contents
	chomp(@Lines);
	close $Map;

	$count=0;
	if(-e $GoaFile){#if the goa file exists
		open ($IPfile, "< $GoaFile")
			or die "Couldn't open input file: $!\n";#open goa file
		@GoaLines=<$IPfile>;#get contents
		#Create hash
		while ($count<@GoaLines) {#for each line of the goa file
			local @Terms=split(/\t/, $GoaLines[$count]);
			if(defined $Annotation{$Terms[$IDColumngoa]}){#if this ID is in the hash
				$Annotation{$Terms[$IDColumngoa]}="$Annotation{$Terms[$IDColumngoa]} $Terms[$GOColumn] $Terms[$ECodeColumn]"; #hash the goa line to the ID
			}
			else{
				$Annotation{$Terms[$IDColumngoa]}="";#clear the string
				$Annotation{$Terms[$IDColumngoa]}="$Terms[$GOColumn] $Terms[$ECodeColumn]";
			}
			
			$count++;
		} #close while loop
		close $IPfile;
	}
	@GoaLines
	
	$count=5;#start at 6th line
	while($count<@Lines){#open while loop
		local @Terms=split(/\t/, $Lines[$count]);
		local $ID=$Terms[8];
		local $Info=$Annotation{$ID};#find the ID
		$OtherFunc=$Terms[-2];#get the function
		$OtherFunc=~s/\|+|\-+/ /g; #replace | or - with a space
		$OtherFunc=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
		$OtherFunc=~ tr/A-Z/a-z/; #convert everything to lower case
		$OtherFunc=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
		$OtherFunc=~ s/^\s+//gm; #remove leading whitespace
		$OtherFunc=~ s/\s+$//; #remove trailing whitespace
		$Start=$Terms[5];#get the Start                          PARSE OTHER FUNCTION
		$Length=$Terms[6];#get the Length
		$Direction=$Terms[7];#get the Direction
		if($Direction eq "F"){#if forward
			$Stop=$Start+$Length;
		}
		else{#else reverse
			$Stop=$Start;
			if($Length<0){
				$Length=$Length*-1;#some CP tables provide a -Length when reverse
			}
			$Start=$Stop+$Length;
		}
		if(defined($Info)){#if its found
			print "$ID\t$Start\t$Stop\t$Info $OtherFunc\n";
		}

		else{#else it was not found
			print "$ID\t$Start\t$Stop\t$OtherFunc\n";
		}
		$count++;
	}#close while loop
	@Lines=();#clear the array
}#close if goa
	





#close $Annotes;
