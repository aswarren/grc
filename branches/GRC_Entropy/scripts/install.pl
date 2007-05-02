#!/usr/bin/perl -w
use Getopt::Std;
#the following is a perl script to compile and copy various components of GRC

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


print "compiling grc_overlap\n";
chdir("../grc_overlap/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_overlap ../../");


if ($status != 0){
	die "grc_overlap did not compile successfully";
}


print "compiling grc_translate\n";
chdir("../grc_translate/");
$status = system("make clean");
$status = $status + system("make");
$status = $status + system("cp ./grc_translate ../../translate/");


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