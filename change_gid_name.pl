#!/bin/perl -w
use strict;

my ($gtf) = @ARGV;
open(IN,"$gtf");
my $usage = "$0 <merged_gtf_file> 
make gid name consistent with transcript id
";
die $usage if @ARGV<1;

my $name_change = $gtf.".use";
open(OUT, "> $name_change") || die "can not open the $name_change\n";
my $n=0;

while(<IN>){
    chomp;
#    if (/^chr[\w]+\t[\w\_]\t[\w\_]+\t[0-9]+\t[0-9]+\t[\w\.]+\t[+-]\t[\w\.]+\t/){
    if (/^(chr[\w]+\t[\w\_]+\t[\w\_]+\t\w+\t\w+\t[\w\.]+\t[+-]\t[\w\.]+)\t.*transcript\_id (\"[\w\_]+\"\;)/){
	$n++;
	print OUT "$1\tgene\_id $2 transcript_id $2\n";

    }
    else {die ("do not match at line $n\n");}
}

