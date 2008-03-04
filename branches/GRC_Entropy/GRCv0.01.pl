#!/usr/bin/perl -w
use Getopt::Std;
use File::Basename;
use Cwd;
use Cwd 'abs_path';
use Time::HiRes qw(gettimeofday);

#this routine retrieves the absolute directory of a file (does not translate links)
sub get_dir{
	my @parms = @_;
	foreach  $p (@parms) {
		$p=abs_path(dirname(glob($p)));
	}
    # Check whether we were called in list context.
    return wantarray ? @parms : $parms[0];
}

sub get_abspath{
        my @parms = @_;
        foreach  $p (@parms) {
                $p=abs_path(dirname(glob($p)))."/".basename($p);
        }
    # Check whether we were called in list context.
    return wantarray ? @parms : $parms[0];
}

sub cat_files {
	my($out_file, @array) = @_;
	open ($out_handle, "> ".$out_file) or die "Couldn't open output file for writing: $!\n";
	select $out_handle;
	foreach $file (@array) {  # indices of @array's elements
		open ($in_handle, "< $file") or die "Couldn't open input file: $!\n";#open file
		while(my $line = <$in_handle>){
			print $line;
		}
	}
	select STDOUT;
}

$start = gettimeofday( );
#Example for running script
# ./GRCv0.01.pl -g AE008687.fna -d testDB
print "\nGRCv0.01 running on @ARGV\n";
my $MinLength=300; #Minimum Length of orfs generated by grc_orfs.  Recommended: >=99
my $Matrix="BLOSUM62"; #Available: BLOSUM80 BLOSUM62 BLOSUM45 PAM70 PAM30
my $NumBHits=10;
my $OntFile="gene_ontology.obo";
my $DBSize=2879860; #effecive size of the DB to adjust e-values accordingly  Set To Zero to use real size


my $CDir=getcwd;#get current working directory
my $BinDir=get_dir($0);#Get the path for current script
#@BinTerms=split(/\//, $BinDir); #split the path for script
#$BinName=$BinTerms[-1]; #set the bin name to be the script name
#$BinDir=~ s/$BinName//g; #remove text
#$BinDir=~ s/\/$//;#remove trailing /
#chomp($BinDir);
#chdir("$BinDir");
#$BinDir=getcwd;#get absolute path for script


$OntFile=$BinDir."/GO/$OntFile";#set absolute ontology filename
my $orfsout ="grc_orfs.out";#variable specifying grc_orf results
my $transdir=$BinDir."/translate/";
my $blastdir=$BinDir."/fsablast/";
my $tempdir=$BinDir."/temp/";
my $transeqout="translate.out";
my $BHName="$tempdir"."bestblast.out";
my $ReferenceName="RefParsed.txt";
my $delim="---------------------------------------------------";
my $DBFile;
my $GFile;
my $GDir;
my $RDir;
my $DBDir;
my $UseGO=0;#boolean variable that tells whether GO is to be used
my $MergeDB=$BinDir."/DB/AutoMerge.faa";
my $TransFile="GCode.txt";#file used for translating
my @AminoFiles;
my %AnnotHash; #annotation complement files stored in hash
my $MergeCommand ="$BinDir"."/scripts/mergeseqannot.pl";
my $ResultDir= $BinDir."/results/";

#chdir("$CDir");

#set MaxBlast FileName (maxblast files are used in grc_overlap calculations)
my $MaxFile="$blastdir"."Max".substr($Matrix,0,1).substr($Matrix,-2).".txt";




#Use command line parameters -g genome.fasta and -d database.faa
getopt('gdrkmh');# get and assign the command line parameters $opt_g $opt_d


unless (-e $opt_g && -e $opt_d) { #check for command line parameter existence
	die "Either $opt_g or $opt_d does not exist\n";
}

if(defined $opt_m){
	$MinLength=$opt_m;
}

if(defined $opt_h){
	$NumBHits=$opt_h;
}


$GDir=get_dir("$opt_g");
$GName=basename("$opt_g");
#@GTerms=split(/\//, $GDir); #split the path for genome
#$GName=$GTerms[-1]; #set the genome name
#$GDir=~ s/$GName//g; #remove text
#$GDir=~ s/\/$//;#remove trailing /
#chomp($GDir);
#chdir("$GDir");
#$GDir=getcwd;
#chdir("$CDir"); #change back to orig wd
$GFile=get_abspath("$opt_g"); #set absolute file name

$GenomeName=(split(/\./, $GName))[0];
$GenomeName=$GenomeName."Min".$MinLength."BH$NumBHits";


if(defined $opt_k){#if the user desires to keep the blast and reference files
	$BHName="$tempdir"."$GenomeName".".bh";
	$ReferenceName=$GenomeName.".ref";
	$opt_k="";#hack shutup perl warning
}



$DBFile=get_abspath("$opt_d");
if (-d $opt_d){#if option d is a directory then merge all sequence and annotation files therein
	chdir("$opt_d");
	$DBDir=abs_path(glob("$opt_d"));#can't submit to get_dir because it will retrieve the directory above
	opendir Direc, "./";
	@contents= readdir Direc; #get the contents of the current directory
	closedir Direc;
	@contentsort = sort @contents;#sort the names
	if(-e "$MergeDB"){
		$status=system("rm -f $MergeDB");
		if($status !=0){
			die "could not remove $MergeDB in danger of dupicating DB\n";
		} 
	}	
	foreach $FileName (@contentsort){#for each file in contents
		unless(-d "$FileName"){#if the file is not a directory
			if(-e "$FileName" && ($FileName=~/faa$/i || $FileName=~/fasta$/i)){#if its an faa file
				push(@AminoFiles, $FileName);
				$FirstName=$FileName;#get name of the file
				$FirstName =~s/.faa$|.fasta$//;#remove extension
				if($FileName=~/faa$/i){
					$AnnotFile=$FirstName.".ptt";
				}
				else {#else the file has to be fasta and its complement is *.goa
					$AnnotFile=$FirstName.".goa";
				}
				if (-e "$DBDir"."/$AnnotFile"){
					$AnnotFile=$DBDir."/$AnnotFile";
					unless(-d "$AnnotFile"){#if its not a directory
						$AnnotHash{$FileName}=$AnnotFile;
					}
				}	
			}
		}#close if not directory
	}#close the foreach loop
	
	cat_files($MergeDB, @AminoFiles);
	
	chdir("$CDir");
}#close if directory


else {#else its not a directory  NOTE Must specify a sequence file *.faa or *.fasta
	$DBDir=get_dir("$opt_d");
	#@DBTerms=split(/\//, $DBDir); #split the path for db
	$DBName=basename(abs_path(glob("$opt_d")));
	$DBFile=get_abspath("$opt_d"); #set path of the db file specified
	unless(-e "$DBFile" && ($DBFile=~/faa$/i || $DBFile=~/fasta$/i)){ #if its the right format and it exists
		die "DB file specified is not compatible or does not exist\n";
	} 
	push(@AminoFiles, $DBName);
	cat_files($MergeDB, ($DBFile));

	$FirstName=$DBName;#get name of the file
	$FirstName =~s/.faa$|.fasta$//;#remove extension
	if($DBName=~/faa$/i){
		$AnnotFile=$FirstName.".ptt";
	}
	else {#else the file has to be fasta and its complement is *.goa
		$AnnotFile=$FirstName.".goa";
	}
		
	if (-e "$DBDir"."/$AnnotFile"){
		$AnnotFile=$DBDir."/$AnnotFile";

		unless(-d "$AnnotFile"){#if its not a directory
			$AnnotHash{$DBName}=$AnnotFile;
		}
	}
	chdir("$CDir");
	
}

$DBFile=$MergeDB;#the database file is now the merged database




if(defined $opt_r && -e $opt_r){#if there is a reference file to compare to
	if($opt_r=~/goa$/i){#if it is a goa reference file
		$UseGO=1;
	}
	$RDir=get_dir("$opt_r");
	#@RTerms=split(/\//, $RDir); #split the path for reference file
	$RName=basename(abs_path(glob("$opt_r"))); #set the reference name
	#$RDir=~ s/$RName//g; #remove text
	#$RDir=~ s/\/$//;#remove trailing /
	#chomp($RDir);
	#chdir("$RDir");
	#$RDir=getcwd;
	#chdir("$CDir"); #change back to orig wd
	$RFile="$RDir"."/$RName"; #set absolute file name
	$status=system("$BinDir"."/scripts/parseRef.pl -i $RFile >$tempdir".$ReferenceName);#parse reference file
	if ($status != 0){
		die "could not find/parse reference file";
	}
}

#change directory to glimmer bin
chdir("$BinDir");

#Run grc_orfs 
print "Running grc_orfs:\n";
$status = system("./grc_orfs $GFile $MinLength $tempdir$orfsout");

if ($status != 0){
	die "grc_orfs did not run successfully";
}


#run grc_translate to translate ORFs to AA
#table 11 use bacterial translation table NO OTHER TABLES SUPPORTED YET
print "\nTranslating sequences.\n";
chdir("$transdir");
$status = system("./grc_translate $transdir$TransFile 11 $tempdir$orfsout $tempdir$transeqout");

if ($status != 0){
	die "grc_translate did not run successfully";
}
chdir("$BinDir");


#create database
print "\nFormatting DB:\n";
$status = system("$blastdir"."formatdb $DBFile");
if($status != 0){
	die "formatdb did not run successfully";
}


#Run FSA-BLAST on AA sequences
print "\nBlasting sequences.\n";
chdir("$blastdir");
$status = system("./blast -d $DBFile -e .001 -i $tempdir$transeqout -m 8 -v 1 -b $NumBHits -M $Matrix -z $DBSize >"."$BHName");

if ($status != 0){
	die "blast did not run successfully";
}
chdir("$BinDir");



#BLAST OUTPUT KEY
#Query id	Subject id	% identity	alignment length	mismatches	gap openings	q. start	q. end	s. start	s. end	e-value	bit score

#AFTER MERGE KEY
# Fields: query id	q. start	q. end	subject id1	subject id2	subject id3	description	organism	% identity	alignment length	subject length	mismatches	gap opens	q. align start	q. align end	s. align start	s. align end	evalue	bit score	frac_filtered
chdir("$DBDir");
foreach $FileName (@AminoFiles){#for each file in contents
	$MergeCommand="$MergeCommand"." $FileName";#add filename to merge command
	local $AnnotFile=$AnnotHash{$FileName};
	if(defined($AnnotFile)){#if its an faa file
		$MergeCommand="$MergeCommand"." $AnnotFile";#add filename to merge command
	}
}#close the foreach loop

$MergeCommand="$MergeCommand $BHName $tempdir$transeqout";
$status=system("$MergeCommand");#run merge
if ($status !=0){
	die "unsuccessful merge\n";
}#close if
chdir("$BinDir");



unless (-e "$ResultDir$GenomeName" && -d "$ResultDir$GenomeName") { #check if the directory exists
	$status = system("mkdir $ResultDir$GenomeName");
	if ($status != 0){
		print "Problem creating $GenomeName directory \n";
	}
}# create the directory if it doesn't exist

$ResultDir= "$ResultDir$GenomeName"."/";

$OverlapCommand=$BinDir."/grc_overlap $BHName $GenomeName $GFile $MaxFile $transdir$TransFile $MinLength";
if($UseGO==1){
	$OverlapCommand=$OverlapCommand." $OntFile";
}

print "$OverlapCommand\n";
#Run GRC_overlap c++ function to remove blast results that are not relevant
#grc_overlap <ParsedBlast> <GenomeName> <Genome.fna>  Output is to cout
print "grc_overlap: adjusting and removing orfs.\n";
chdir("$ResultDir");
$status = system($OverlapCommand);
if ($status != 0){
	die "grc_overlap did not run successfully";
}

$elapsed = gettimeofday( ) - $start;
$Minutes=$elapsed/60;
print "Running time: $Minutes minutes\n";
#Run GRC_compare:  for comparing output of the GRC to a reference file provided

if(defined $opt_r){#if there is a reference file to compare to
	print "grc_compare: comparing to reference file $opt_r\n";
	$CompareCommand=$BinDir."/grc_compare $tempdir"."$ReferenceName $GenomeName".".Pos $GenomeName".".Neg ./KnockList.txt $MinLength";
	if($UseGO==1){
		$CompareCommand=$CompareCommand." $OntFile";
	}
	$CompareCommand=$CompareCommand." >$GenomeName".".compare";
	print "$CompareCommand\n";
	$status = system($CompareCommand);
	#$ParseCommand=$BinDir."/scripts/ParseCompare.pl ./$GenomeName".".compare >$GenomeName".".comparsed";
	#system($ParseCommand);
	if ($status != 0){
		die "grc_compare did not run successfully";
	}
}
	
	

#move other results files
#$status= system("mv ./KnockList.txt $ResultDir");#move the knocklist to the compare directory
#$status = $status + system("mv ./$Negatives $ResultDir$NegParsed");#move the negatives file
#$status = $status + system("mv ./$Positives $ResultDir$PosParsed");#move the positives file
#if ($status != 0){
#		print "Problem moving additional $GenomeName results \n";
#	}




