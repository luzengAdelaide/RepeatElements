#!/usr/bin/perl -w
use strict;
 
my ($file1,$file2) = @ARGV;
my $usage = "usage: $0 <longest genes> <merged transcripts> 
Extract the gene detail (exon and cds) from longest genes output.
";

open(INA,"$file1");
open(INB,"$file2");
my $out = "single.transcript";
open(OUT,">$out");
my ($go);
my %hash;

while (<INA>) {
    chomp;
    my @data = split("\t",$_);
    $go = $data[9];
    $hash{$go} = 1;
}

while(<INB>) {
    chomp;
    my @tmp = split("\t",$_);
    print OUT $_,"\n" if exists $hash{$tmp[9]};
}

close INA;
close INB;
