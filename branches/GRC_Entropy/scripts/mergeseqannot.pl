#!/usr/bin/perl -w
use Getopt::Std;
#the following is a perl script to merge AA sequence and annotation files 
#These must be provided in (annotation, sequence) pairs e.g. "./mergeseqannot.pl a1.faa a1.ptt a2.faa a2.ptt" with the last parameter being the output file
#e.g. "./mergeseqannot.pl a1.fasta a1.goa a2.fasta a2.goa"

#BLAST OUTPUT KEY
#Query id	Subject id	% identity	alignment length	mismatches	gap openings	q. start	q. end	s. start	s. end	e-value	bit score

#AFTER MERGE KEY
# Fields: query id	q. start	q. end	subject id1	subject id2	subject id3	description	organism	% identity	alignment length	subject length	mismatches	gap opens	q. align start	q. align end	s. align start	s. align end	evalue	bit score	frac_filtered

sub get_position_info{
	local($source) =@_;
	local @SubTerms=split(/\_/, $source);
	$SubTerms[0]=~ s/^>//gm; #remove leading carrot if it exists
	return ($SubTerms[0], $SubTerms[1], $SubTerms[2]);
}

sub get_replicon_info {
	local($source) = @_;
	local $ReturnID="";
	local $ReturnInfo="";
	if($source=~/\|/){
		@IDSplit2=split(/\|/,$source);
		$ReturnID=$IDSplit2[0];
		$ReturnInfo=$source;
		$ReturnInfo=~s/$LookupID//;#remove ID with ID_start_stop
	}
	else{
		$ReturnID=$source;
		$ReturnInfo="";
	}
	return ($ReturnID, $ReturnInfo);
}

sub get_subject_id{
	local($source) =@_;
	if($source=~/\|/){#if the ID contains a seperator
		$source=~ s/^\|+//gm; #remove leading seperator '|'
		@IDSplit=split(/\|/,$source);
		if($IDSplit[0] eq "gi"){
			$source=$IDSplit[1];	
		}
		else{
			$source=$IDSplit[0];
		}
	}
	return $source;
}

my %Annotation =(); #hash to store annotation information
my $NumParam=@ARGV;#get the number of parameters
$NumParam=$NumParam-2;#subtract off output parameter
$TotalCount=0;
print "Merging sequence and annotation files\n";
my $blast_output=$ARGV[-2];
my $orf_file=$ARGV[-1];

#rename blast output
$status=system("mv -f $blast_output $blast_output".".orig");
if($status !=0){
	die "could not copy file in merge procedure. in danger of losing blast output\n";
}

$blast_input="$blast_output".".orig";

#read in blast input
open ($blastfile, "< ".$blast_input) or die "Couldn't open blast output for parsing: $!\n";
my @BlastLines=<$blastfile>;#get contents
chomp(@BlastLines);

#create hash table based on query id for blast results
my %BlastHash;
$count=0;
while ($count < scalar @BlastLines){
	local @Terms=split(/\t/, $BlastLines[$count]);
	local $TempID=get_subject_id($Terms[1]);

	$BlastHash{$TempID}="INITIAL";
	$count=$count+1;
}



$a=0;
while($a<$NumParam){ 

	if($a+1<$NumParam && $ARGV[$a+1]=~/ptt$/i){#open parameter loop
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
		
		

		close $IPfile;
		@Lines=();
	}

	if($a<$NumParam && $ARGV[$a]=~/faa$/i){
		$RecordCount=0;
		$InHandle=$ARGV[$a];#set the input file name
		open ($IPfile, "< $InHandle")
			or die "Couldn't open input file: $!\n";#open faa file

		local $PreviousLine;
		local $CurrentLength=0;
		local $EndLine;
		local $PreviousID="!NOBLASTRESULT!";
		while(<$IPfile>){#open faa for loop
			if($_=~/^>/){#if aa record header
				if($CurrentLength!=0 && $PreviousID ne "!NOBLASTRESULT!"){
					$BlastHash{$PreviousID}="$PreviousLine"."$CurrentLength";
					$PreviousID="!NOBLASTRESULT!";
				}
				$CurrentLength=0;
				$RecordCount++;
				local @Terms=split(/\|+/, $_);#break up line to get ID
				local $SubjectID= $Terms[1];# in faa file format expected is ">gi|number|etc."
	
				local $BLine=$BlastHash{$SubjectID};#look up blast results
				if(defined($BLine) && $BLine eq "INITIAL"){
					local $Info=$Annotation{$SubjectID};#use 2nd element (should be ID) to get ptt line
					local $GeneName = "-"; #the gene name from the ptt file
					local $Synonym ="-"; #the synonym code from the ptt file
					if(defined($Info)){#if its found in the ptt hash
						local @AnnTerms=split(/\t/,$Info);#split by tabs
						chomp @AnnTerms;
						$GeneName =$AnnTerms[4]; #the gene name from the ptt file
						$Synonym =$AnnTerms[5]; #the synonym code from the ptt file
						$Function="$AnnTerms[-1]";#assign function description to the product and gene name
					}#close if in hash
					else{#else not in hash
						$Function=$Terms[-1];#everything after last | in faa header
					}
					$Function=~s/\|+|\-+/ /g; #replace | or - with a space
					$Function=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
					#$Function=~ tr/A-Z/a-z/; #convert everything to lower case
					#$Function=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
					$Function=~ s/^\s+//gm; #remove leading whitespace
					$Function=~ s/\s+$//; #remove trailing whitespace

					$PreviousID=$SubjectID;
					#local @BTerms=split(/\t/, $BLine);
					#local @SubTerms=split(/\_/, $BTerms[0]);
					#               subject id1   subject id2   subject id3   description   organism
					$PreviousLine="$SubjectID\t$GeneName\t$Synonym\t$Function\t$InHandle\t";
				}
			}#close if header
			else{#if its not a header
				chomp;
				$CurrentLength+=length($_);
				#print $OPfile $_."\n";
			}
		}#close while faa open
		if($CurrentLength!=0 && $PreviousID ne "!NOBLASTRESULT!"){
			$BlastHash{$PreviousID}= "$PreviousLine"."$CurrentLength";
			$PreviousID="!NOBLASTRESULT!";
		}
		%Annotation=();#clear hash
		print "$InHandle: $RecordCount records\n";
		$TotalCount=$RecordCount+$TotalCount;
	}

		#else{#else annotation sequence file mismatch
		#	die "error combining files\n";
		#}
		
#BLAST OUTPUT KEY
#Query id	Subject id	% identity	alignment length	mismatches	gap openings	q. start	q. end	s. start	s. end	e-value	bit score

#AFTER MERGE KEY
# Fields: query id	q. start	q. end	subject id1	subject id2	subject id3	description	organism	% identity	alignment length	subject length	mismatches	gap opens	q. align start	q. align end	s. align start	s. align end	evalue	bit score	frac_filtered





	if($a+1<$NumParam && $ARGV[$a+1]=~/goa$/i){
		$IDPosition=1;
		$GOPosition=4;
		$ECodePosition=6;


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
		close $IPfile;
		@Lines=();
	}

	if($a<$NumParam && $ARGV[$a]=~/fasta$/i){
		local $PreviousLine;
		local $CurrentLength=0;
		local $EndLine;
		local $PreviousID="!NOBLASTRESULT!";
		$RecordCount=0;
		$InHandle=$ARGV[$a];#set the input file name
		open ($IPfile, "< $InHandle")
			or die "Couldn't open input file: $!\n";#open fasta file
		while(<$IPfile>){#open faa for loop
			chomp;
			if($_=~/^>/){#if aa record header
				if($CurrentLength!=0 && $PreviousID ne "!NOBLASTRESULT!"){
					$BlastHash{$PreviousID}="$PreviousLine"."$CurrentLength";
					$PreviousID="!NOBLASTRESULT!";
				}
				$CurrentLength=0;
				$RecordCount++;
				local @Terms=split(/\s+/, $_);#break up line on whitespace to get ID
				local $ID=$Terms[0];
				$ID=~s/>//;#remove '>'
				if($ID=~/\|/){#if the ID contains a seperator
					@IDSplit=split(/\|/,$ID);
					$ID=$IDSplit[0];
					#$ID=~s/\|$SecondID//g;#remove everything after seperator
					#print "$ID\n";
				}
				local $BLine=$BlastHash{$ID};#look up blast results
				if(defined($BLine)&& $BLine eq "INITIAL"){
					$Function=$_;
					$Function=~s/>$ID//; #remove ID from the rest of the annotation
					$Function=~s/\|+|\-+/ /g; #replace | or - with a space
					$Function=~s/\(+|\)+|\[+|\]+|\{+|\}+|,+|\/+|\\+|\,+|\-+|\'+|\;+|\t+/ /g; #replace brackets and punct with space
					#$Function=~ tr/A-Z/a-z/; #convert everything to lower case
					#$Function=~s/ protein +| and +| the +/ /g; #remove undesirable annotations
					$Function=~ s/^\s+//gm; #remove leading whitespace
					$Function=~ s/\s+$//; #remove trailing whitespace
					local $Info=$Annotation{$ID};#use ID to get GO annotations
					local $GeneName = "-"; #the gene name 
					local $Synonym ="-"; #the synonym code 

					if(defined($Info)){#if its found
						#local @AnnTerms=split(/\t+/,$Info);#split by tabs
						#chomp @AnnTerms;
						#print $OPfile ">".$ID."\t$Info $Function\t$InHandle\n";#print ">ID ProductDescription FileName"
						$Function="$Info $Function";
					}#close if in hash

					$PreviousID=$ID;

					#           query id   q. start   q. end   subject id1   subject id2   subject id3   description   organism   % identity   alignment length
					$PreviousLine="$ID\t$GeneName\t$Synonym\t$Function\t$InHandle\t";
				}
				
			}#close if header
			else{#if its not a header
				$CurrentLength+=length($_);
				#print $OPfile $_."\n";
			}
		}#close while fasta open
		if($CurrentLength!=0 && $PreviousID ne "!NOBLASTRESULT!"){
			$BlastHash{$PreviousID}="$PreviousLine"."$CurrentLength";
			$PreviousID="!NOBLASTRESULT!";
		}
		%Annotation=();#clear hash
		print "$InHandle: $RecordCount records\t";
		$TotalCount=$RecordCount+$TotalCount;
	}#close if fasta file
	$a++; #increment counter
}#close while loop

if($TotalCount==0){#else annotation sequence file mismatch
	die "error in merging annotation database files in mergeseqannot.pl. No files for database.\n";
}


#open write to old blast filename
open ($op_handle, "> ".$blast_output) or die "Couldn't open output file: $!\n";
open ($orf_handle, "< ".$orf_file) or die "Couldn't open orfs file: $!\n";

$count=0;
my @BTerms=split(/\t/, $BlastLines[$count]);
my ($BlastID, $BlastReplicon)=get_replicon_info($BTerms[0]);


while (<$orf_handle>){
	$Line=$_;
	chomp($Line);
	if($Line=~/^>/){#if aa record header
		$Line=~ s/^>//gm; #remove leading carrot if it exists

		@OrfTerms=split(/\t/, $Line);
		($LookupID, $RepliconInfo)=get_replicon_info($OrfTerms[0]);
		
		
	
		if($LookupID ne $BlastID){#the orf has no hit
			local ($QueryID, $Start, $Stop) = get_position_info($LookupID);
			$QueryID="$QueryID"."$RepliconInfo";
			print $op_handle "$QueryID\t$Start\t$Stop\t"."No_hits\n";
		}

		else{
			while($BlastID eq $LookupID){
				local ($QueryID, $Start, $Stop) = get_position_info($BlastID);
				$QueryID="$QueryID"."$BlastReplicon";
				local $SubjectID=get_subject_id($BTerms[1]);
	
				local $SubjectInfo=$BlastHash{$SubjectID};
				if(defined($SubjectInfo)){
					($SID, $GeneName, $Synonym, $Function, $Organism, $SubjectLength)=split(/\t/,$SubjectInfo);
				}
				else {
					print "$SubjectID\n";
					$SID=$GeneName=$Synonym=$Function=$Organism="-";
					$SubjectLength=0;
				}
				#                query id   q. start   q. end   subject id1   subject id2   subject id3   description   organism   % identity   alignment length
				$PrintLine="$QueryID\t$Start\t$Stop\t$SubjectID\t$GeneName\t$Synonym\t$Function\t$Organism\t$BTerms[2]\t$BTerms[3]\t";
				#          mismatches   gap opens   q. align start   q. align end   s. align start   s. align end   evalue   bit score
				$PrintLine="$PrintLine"."$SubjectLength\t$BTerms[4]\t$BTerms[5]\t$BTerms[6]\t$BTerms[7]\t$BTerms[8]\t$BTerms[9]\t$BTerms[10]\t$BTerms[11]\n";
				print $op_handle "$PrintLine";
				$count=$count+1;
				if($count < scalar @BlastLines){
					@BTerms=split(/\t/, $BlastLines[$count]);
					@Results=get_replicon_info($BTerms[0]);
					$BlastID=$Results[0];
					$BlastReplicon=$Results[1];
				}
				else{
					$BlastID="!!NONE!!";
				}
			}#close while
		}
	}
}




print "TotalCount: $TotalCount \n\n";

	
close $op_handle;





