#!/usr/bin/perl 

my ($file) = @ARGV;
my $usage = "usage: $0 <merged_transcript> 
Acquire the longest genes from merged transcripts.
";

open(IN,"$file");
my $out = "longest.".$file;
open(OUT,">$out");
my %hash;
my $chr;
my @tmp;
my $info;
my $val;
my $i=0;
my @gene_id;
my @chr;
my $test;
my @array1="";
my $chr_pre="";
my $chr_aft="";
my @all_keys;

while(<IN>) {
    @tmp = split("\t",$_);
    $chr=$tmp[0].$tmp[8]; 
    next if ($tmp[0] =~ /\_/ );
    $hash{$chr} = $hash{$chr}.$_;
    $chr_aft = $chr;
    if($chr_aft ne $chr_pre){
	push (@all_keys,$chr_aft);
	$chr_pre=$chr_aft;
    }else{
	$chr_pre=$chr_aft;
    }
}

foreach $val  (@all_keys) {
    my ($min,$max)=(999999999,0);
    @gene_id = split("\n",$hash{$val});
    foreach $test_tmp (@gene_id) {
	@test = split("\t",$test_tmp);
	$min=$test[3] if($test[3]<$min);
	$max=$test[4] if($test[4]>$max);
    }
    print OUT "$test[0]\t$test[1]\t$test[2]\t$min\t$max\t$test[5]\t$test[6]\t$test[7]\t$test[8]\t$test[9]\n";
}
