#!/usr/bin/perl -w
use strict;

# please change the column number of @head and @tail based on the overlapped datasets

my ($filename)=@ARGV;
my $usage = "usage: $0 <overlapped region> > <calculated result>
";

open(IN,"$filename");
my (@head,@tail,@length);

while(<IN>){
    chomp;
    my $n++;
    my @ref=split("\t",$_);
    push @head,$ref[1];
    push @tail,$ref[2];
} 
close IN;

my $length=0;
for(my $i=0;$i<@head;$i++){
     $length ++;
     $length[$i]=$tail[$i]-$head[$i];
#     print "the length is $length\n";
}


my $total=0;
for(my $i=0;$i<@head;$i++){
    $total += $length[$i]+1;
}
print "the total number is $total\n";
