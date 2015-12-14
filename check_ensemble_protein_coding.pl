#!/usr/bin/perl -w 


## written by selfless sporting steady gentle man, i.e. df_cao.

use strict;
my $usage = "usage: ./check_ensemble_protein_coding.pl <genome.protein_coding>\n
Delete bad gene sequence lines in genome with gtf format
";
die $usage if @ARGV<1;

my ($genome) = @ARGV;
open(GENOME, $genome) || die "can not open the file $genome\n";
my $right_genome = $genome . ".check";
open(OUT, "> $right_genome") || die "can not open the $right_genome\n";

my (@gencode,@temp_gene);
my ($i,$id_last,$flag, $flag_1,$have_exon,$have_CDS)=(0,0,0,0,0,0);

print "waiting!!!processing genome.\n";
while(<GENOME>)
{
    chomp;
    $gencode[$i++] = $_;
}
$i=1;
foreach(@gencode)
{
    chomp;
    my $l_now = $_;
    $flag_1=0;
    if ($l_now =~ /chr[\w]+\_[\w]+/) {
	next;
    }
    if($l_now =~ /.*transcript_id\s+\"(\w+\.?(\w+)?\-?(\w+)?)\"/) 
    {  
	my $id_now =$1;
	if($flag == 0)
	{
	    $id_last = $id_now;
	    $temp_gene[0] = $l_now;
	    $flag=1;
	    next;
	}
	#store a gene in @temp_gene;
	if($id_now eq $id_last){$temp_gene[$i++] = $l_now;next;}
	elsif($id_now ne $id_last)
	{
	    unshift(@gencode, $l_now); #notice:@gencode has been changed!!!;
	    for(my $k=0;($k<=$#temp_gene && $flag_1==0);++$k)
	    {
		#if($id_now eq "ENST00000569434"){print "temp_gene = $#temp_gene,[$k]=$temp_gene[$k]\n";}
		if(($temp_gene[$k])&&($temp_gene[$k] =~ /(start|stop)_codon\s+(\d+)\s+(\d+)/))
		{
		    my ($promoter,$st_now, $en_now)=($1,$2,$3);
		    if(abs($st_now-$en_now) != 2){&clear_temp_gene;last;}
		    else
		    {
			my $flag_2=0;
			if($promoter eq "start")
			{
			    for(my $j=$k+1;$j<=$#temp_gene;++$j)
			    {	
				if($temp_gene[$j] =~ /start_codon\s+\d+\s+\d+/){&clear_temp_gene;last;}
				#elsif(($temp_gene[$j] =~ /(CDS|exon)\s+(\d+)\s+(\d+)/)&&($flag_2==0)){$flag_2=1;}
				elsif(($temp_gene[$j] =~ /stop_codon\s+\d+\s+\d+/)&&($flag_2==0)){$flag_2=1;}
				elsif(($temp_gene[$j] =~ /stop_codon\s+\d+\s+\d+/)&&($flag_2==1)){&clear_temp_gene;last;}
				elsif($temp_gene[$j] =~ /exon\s+\d+\s+\d+/){$have_exon=1;}
				elsif($temp_gene[$j] =~ /CDS\s+\d+\s+\d+/){$have_CDS=1;}
				
#if(($j == $#temp_gene)&&($flag_2 == 1)&&($have_CDS == 1)&&($have_exon == 1))
				if(($j == $#temp_gene)&&($flag_2 == 1)&&($have_exon == 1))
				{
				    for(my $q=0;$q<=$#temp_gene;++$q)
				    {print OUT "$temp_gene[$q]\n";}
				    #print ".";
				    &clear_temp_gene;last;
				}
			    }
			    last;
			}
			elsif($promoter eq "stop")
			{
			    #print "$l_now\n";
			    for(my $j=$k+1;$j<=$#temp_gene;++$j)
			    {
				if($temp_gene[$j] =~ /stop_codon\s+\d+\s+\d+/){&clear_temp_gene;last;}
				#elsif(($temp_gene[$j] =~ /(CDS|exon)\s+(\d+)\s+(\d+)/)&&($flag_2==0)){$flag_2=1;}
				elsif(($temp_gene[$j] =~ /start_codon\s+\d+\s+\d+/)&&($flag_2==0)){$flag_2=1;}
				elsif(($temp_gene[$j] =~ /start_codon\s+\d+\s+\d+/)&&($flag_2==1)){&clear_temp_gene;last;}
				elsif($temp_gene[$j] =~ /exon\s+\d+\s+\d+/){$have_exon=1;}
				elsif($temp_gene[$j] =~ /CDS\s+\d+\s+\d+/){$have_CDS=1;}
				
#if(($j == $#temp_gene)&&($flag_2 == 1)&&($have_CDS == 1)&&($have_exon == 1))
				if(($j == $#temp_gene)&&($flag_2 == 1)&&($have_exon == 1))
				{
				    for(my $q=0;$q<=$#temp_gene;++$q)
				    {print OUT "$temp_gene[$q]\n";}
				    #print ".";
				    &clear_temp_gene;last;
				}
			    }
			    last;
			}
		    } #if(abs($st_now-$en_now) == 2);
		}#if($temp_gene[$k] =~ /(start|stop)_codon\s+(\d+)\s+(\d+)/);
		elsif($temp_gene[$k] =~ /exon\s+(\d+)\s+(\d+)/){$have_exon=1;next;}
		elsif($temp_gene[$k] =~ /CDS\s+(\d+)\s+(\d+)/){$have_CDS=1;next;}
		else{&clear_temp_gene;last;}
	    }
	    &clear_temp_gene;
	}
	$id_last = $id_now;
    }
}

#clear @temp_gene, $flag_1=1;
sub clear_temp_gene
{
    @temp_gene = ( );
    $flag=0;$i=1;$flag_1=1;
    ($have_exon,$have_CDS)=(0,0);
}
