#!/usr/bin/perl -w
use strict;


my $dir = "Dec_dis/";
opendir(DIR,$dir) or die $!;
#@docs = grep (/[\w\_]+/, readdir(DIR));

#foreach $dir (@docs){ 
while (my $file=readdir(DIR)) {
    open (IN,$file) or die "could not open $file\n";
    my $out=$file."_dec";
    open (OUT,">$out");
    my (@chr,@start,@stop,@feature,@strand,@gid);
    my $i=0;
    
    while(<IN>){
	chomp;
	my @data = split("\t",$_);
	push @chr, $data[0];
	push @start, $data[1];
	push @stop, $data[2];
	push @feature, $data[3];
	push @strand, $data[4];
	push @gid, $data[5];
    }
    for($i=0;$i<@chr;$i++){
	if (($stop[$i]-$start[$i]+1) >=25) {
	    print OUT "$chr[$i]\t$start[$i]\t$stop[$i]\t$feature[$i]\t$strand[$i]\t$gid[$i]\n";
	}
	else {
	    next;
	}
    }
}
