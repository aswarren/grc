#!/usr/bin/perl -w
use Getopt::Std;
use File::Basename;
use Cwd;
#the following is a perl script to run annotations in a leave one out fashion
#this script assumes that all files for an organism are named the same eg. a.faa a.ptt a.goa a.fasta. a.CP a.fna
#this script also assumes you are running it from the scripts folder in GRC and specify the genome, reference and db folder relative to that
#../DB/;../genomes/;../reference/;../DB/
$MinSize=99;
$BHits=10;

my $MainScript = "../GRCvNo_Blast.pl";
my $CDir=getcwd;#get current working directory
getopt('gdrmbs');# get and assign the command line parameters $opt_g $opt_d

unless (-e $opt_g && -e $opt_d) { #check for command line parameter existence
	die "Either $opt_g or $opt_d does not exist\n";
}
my $GenomeFolder=$opt_g;
my $DBFolder=$opt_d;
my $ReferenceFolder="none";
if(defined $opt_r){
	$ReferenceFolder=$opt_r;
	chdir("$opt_r");
	$RDir=getcwd;
	opendir Direc, "./";
	@rcontents= readdir Direc; #get the contents of the current directory
	closedir Direc;
	@rcontentsort = sort @rcontents;#sort the names
	chdir("$CDir"); #change back to orig wd
}

if(defined $opt_s){
	$MainScript=$opt_s;
}

if(defined $opt_m){
    $MinSize=$opt_m;
}

if(defined $opt_b){
    $BHits=$opt_b;
}

chdir("$opt_d");
$DBDir=getcwd;
opendir Direc, "./";
@dcontents= readdir Direc; #get the contents of the current directory
closedir Direc;
@dcontentsort = sort @dcontents;#sort the names
chdir("$CDir"); #change back to orig wd

chdir("$opt_g");
$GDir=getcwd;
opendir Direc, "./";
@gcontents= readdir Direc; #get the contents of the current directory
closedir Direc;
@gcontentsort = sort @gcontents;#sort the names
chdir("$CDir"); #change back to orig wd
my $RefFileName;

foreach $Genome(@gcontentsort){#for each file in the genomes folder
	if($Genome=~/.fna/){#if its an fna file
		$GenomeName=$Genome;
		$GenomeName=~s/.fna//;#remove extension
		chdir("$opt_d");#change to the DB file
		foreach $DFile(@dcontentsort){#leave out the current genome to be annotated
			if($DFile=~/$GenomeName/){
				if($DFile=~/.goa/){
					$status=system("mv ./$DFile ./$GenomeName".".g");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if goa
				if($DFile=~/.fasta/){
					$status=system("mv ./$DFile ./$GenomeName".".fe");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if fasta
				if($DFile=~/.faa/){
					$status=system("mv ./$DFile ./$GenomeName".".fn");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if faa
				if($DFile=~/.ptt/){
					$status=system("mv ./$DFile ./$GenomeName".".p");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if ptt
			}
		}
		chdir("$CDir"); #change back to orig wd

		foreach$RFile(@rcontentsort){
			if($RFile=~/$GenomeName/){
				if($RFile=~/.goa/){
					$RefFileName=$RFile;
				}
				elsif($RFile=~/.ptt/){
					$RefFileName=$RFile;
				}
			}
		}

		$RunCommand="$MainScript -g $GenomeFolder"."/$Genome -d $DBFolder -m $MinSize -h $BHits -k $GenomeName";
		if(defined $opt_r){
			$RunCommand=$RunCommand." -r $ReferenceFolder"."/$RefFileName >$GenomeName$MinSize$BHits".".TimeSize";
		}

		print "Leave-one-out running $RunCommand\n";
		$status=system("$RunCommand");
		if($status !=0){
				die "failed on $RunCommand\n";
			}
		
		chdir("$opt_d");#change to the DB file

		foreach $DFile(@dcontentsort){#change the filenames back
			if($DFile=~/$GenomeName/){
				if($DFile=~/.goa/){
					$status=system("mv ./$GenomeName".".g ./$DFile" );
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if goa
				if($DFile=~/.fasta/){
					$status=system("mv ./$GenomeName".".fe ./$DFile ");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if fasta
				if($DFile=~/.faa/){
					$status=system("mv ./$GenomeName".".fn ./$DFile");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if faa
				if($DFile=~/.ptt/){
					$status=system("mv ./$GenomeName".".p ./$DFile");
					if($status !=0){
						die "could not rename $DFile\n";
					}
				}#close if ptt
			}
		}
		chdir("$CDir"); #change back to orig wd
	}#close if fna
}#close for each genome

		




