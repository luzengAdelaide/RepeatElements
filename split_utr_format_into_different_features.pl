#!/usr/bin/perl -w 
##--------------------------------------------------------------------------
## Description:
## split UTR format protein coding into different features.
##--------------------------------------------------------------------------
## written by selfless sporting steady gentle man, i.e. df_cao. 
##*************************************************************************** 
## Revervion:0.2 ## 2013/8/16 ## by df_cao ## applied to <ensemble.gtf>.
##***************************************************************************
use strict;

my $usage = "usage: perl split_utr_format_into_different_features.pl <UTR_format> <nonsex_chrs>\n";
die $usage unless @ARGV == 2;

my ($genome,$nonsex_chrs)=@ARGV;
open(GENOME, $genome) || die "cannot open the file $genome\n";

open(U5,"> 5_UTR_exon");
open(U5I,"> 5_UTR_intron");
open(START,"> start_codon");
open(CE,"> CDS_exon");
open(CI,"> CDS_intron");
open(STOP,"> stop_codon");
open(U3,"> 3_UTR_exon");
open(U3I,"> 3_UTR_intron");
open(TRANS, "> protein_coding_transcript");

my($last_st,$last_en)=(0,0);
my ($chr_last,$last_id,$flag) = (0,0,0);
my $gene_st;
my ($max,$min)=(0,0);

#replace position of features in '-' sequence.
my (@gencode,@temp_gene);
my ($i_temp,$flag_temp,$flag_temp_1)=(0,0,0);
my $id_last;

while(<GENOME>)
{
    chomp;
    my $line = $_;
    if($line =~ /^chr(\w+)\s+\w+\s+(\w+)\s+(\d+)\s+(\d+)\s+.*\s+(\-|\+)\s+.*transcript_id\s+\"(\w+\.?\w+)\"/)
    {#print "here\n";
	my ($chr,$type,$st,$en,$op,$id) = ($1,$2,$3,$4,$5,$6);
	if($chr eq 'X'){$chr = ($nonsex_chrs+1);}
	elsif($chr eq 'Y'){$chr = $nonsex_chrs+2;}
	elsif($chr =~ /[a-zA-Z]/){next;}
	
	if($chr ne $chr_last){print "waiting...chr$chr is processing!\n";$flag = 0;}
	if($last_id ne $id)
	{
	    $flag = 0;
	    if($last_en != 0)
	    {print TRANS "chr$chr_last\ttranscript\t$min\t$max\t$op\ttranscript_id\t$id\n";}
	}
	
	if($flag == 0)
	{
	    $gene_st = $st;
	    ($max,$min)=($en,$st);
	    if($type eq "5_UTR"){print U5 "chr$chr\t5_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 1;}
	    elsif($type =~ /start/){print START "chr$chr\tstart\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 2;}
	    elsif($type eq "3_UTR"){print U3 "chr$chr\t3_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 6;}
	    elsif($type =~ /stop/){print STOP "chr$chr\tstop\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 7;}
	}
	elsif($flag == 1)
	{
	    &update($st,$en);
	    if($type eq "5_UTR")
	    {
		print U5 "chr$chr\t5_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print U5I "chr$chr\t5_UTR_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	    elsif($type =~ /start/)
	    {
		#my ($in_st,$in_en) = ($last_en+1,$st-1);
		#print U5I "chr$chr\t5_UTR_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
		print START "chr$chr\tstart\t$st\t$en\t$op\ttranscript_id\t$id\n";
		$flag = 2;
	    }
	}
	elsif($flag == 2)
	{
	    &update($st,$en);
	    if($type =~ /stop/){print STOP "chr$chr\tstop\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 4;}
	    elsif($type =~ /CDS/)
	    {print CE "chr$chr\tCDS_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 3;}
	}	
	elsif($flag == 3)
	{
	    &update($st,$en);
	    if($type =~ /stop/){print STOP "chr$chr\tstop\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 4;}
	    elsif($type =~ /CDS/)
	    {
		print CE "chr$chr\tCDS_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print CI "chr$chr\tCDS_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	}
	elsif($flag == 4)
	{
	    &update($st,$en);
	    if($type eq "3_UTR")
	    {
		#my ($in_st,$in_en) = ($last_en+1,$st-1);
		#print U3I "chr$chr\t3_UTR_intron\t$in_st\t$in_en\n";
		print U3 "chr$chr\t3_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		$flag = 5;
	    }
	}
	elsif($flag == 5)
	{
	    &update($st,$en);
	    if($type eq "3_UTR")
	    {
		print U3 "chr$chr\t3_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print U3I "chr$chr\t3_UTR_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	}
	elsif($flag == 6)
	{
	    &update($st,$en);
	    if($type eq "3_UTR")
	    {
		print U3 "chr$chr\t3_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print U3I "chr$chr\t3_UTR_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	    elsif($type =~ /stop/){print STOP "chr$chr\tstop\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 7;}
	}
	elsif($flag == 7)
	{
	    &update($st,$en);
	    if($type =~ /start/){print START "chr$chr\tstart\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 9;}
	    elsif($type =~ /CDS/){print CE "chr$chr\tCDS_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 8;}
	}	
	elsif($flag == 8)
	{
	    &update($st,$en);
	    if($type =~ /start/){print START "chr$chr\tstart\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 9;}
	    elsif($type =~ /CDS/)
	    {
		print CE "chr$chr\tCDS_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print CI "chr$chr\tCDS_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	}
	elsif($flag == 9)
	{
	    &update($st,$en);
	    if($type eq "5_UTR"){print U5 "chr$chr\t5_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";$flag = 10;}
	}
	elsif($flag == 10)
	{
	    &update($st,$en);
	    if($type eq "5_UTR")
	    {
		print U5 "chr$chr\t5_UTR_exon\t$st\t$en\t$op\ttranscript_id\t$id\n";
		my ($in_st,$in_en) = ($last_en+1,$st-1);
		print U5I "chr$chr\t5_UTR_intron\t$in_st\t$in_en\t$op\ttranscript_id\t$id\n";
	    }
	}
	$last_st = $st;
	$last_en = $en;
	$last_id = $id;
	$chr_last = $chr;
    }
}
#system("rm $temp_utr");

####################-----subroutine-----####################
sub update()
{
    my ($a,$b)=@_;
    if($a > $max){$max = $a;}
    if($b > $max){$max = $b;}
    if($a < $min){$min = $a;}
    if($b < $min){$min = $b;}
}
