#!/usr/bin/perl -w
use Getopt::Std;
use File::Basename;
use Cwd;
use Cwd 'abs_path';
#the following is a perl script to compile and copy various components of GRC

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
my $CDir=getcwd;#get current working directory
my $BinDir=get_dir($0);#Get the path for current script
chdir("$BinDir");

print "compiling grc_orfs\n";
chdir("../src/grc_orfs/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_orfs ../../");


if ($status != 0){
	die "grc_orfs did not compile successfully";
}


print "compiling fsablast\n";
chdir("../fsablast/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./blast ../../fsablast/");
$status = $status + system("cp ./formatdb ../../fsablast/");


if ($status != 0){
	die "blast did not compile successfully";
}


print "compiling grc_annotate\n";
chdir("../grc_annotate/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_annotate ../../");


if ($status != 0){
	die "grc_annotate did not compile successfully";
}


print "compiling grc_translate\n";
chdir("../grc_translate/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_translate ../../");


if ($status != 0){
	die "grc_translate did not compile successfully";
}

print "compiling grc_compare\n";
chdir("../grc_compare/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_compare ../../");


if ($status != 0){
	die "grc_compare did not compile successfully";
}

chdir("$CDir");
