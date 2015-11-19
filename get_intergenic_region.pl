#!/usr/bin/perl -w

##---------------------------------------------------------------------  
## Description:
## add intergenic region information in genomes with gtf format.
##---------------------------------------------------------------------
## written by zhiqiang Hu
 

use strict;

my $usage = "$0 <input_file> <genome_directory>
";
die $usage if @ARGV!=2;
my ($gtf,$gd) = @ARGV;

my @chr = 1..22;
$chr[22] = 'X';
$chr[23] = 'Y';
#temp directories for intermediate computating
system("mkdir tmp");
system("mkdir tmp/annotation");
system("mkdir tmp/prediction");

#cut gtf files into chromosomes
print STDERR "cut gtf files into chromosomes......";
for(my $i=0;$i<@chr;$i++){
    system("grep -P \"^chr$chr[$i]\t\" $gtf >tmp/annotation/$chr[$i]");
}
print STDERR "OK!\n";

for(my $i=0;$i<@chr;$i++){
    print STDERR "chr$chr[$i] processing...";
    my $ann_gtf = "tmp/annotation/".$chr[$i];
    my ($Vec, $num) = @{Read_GTF($ann_gtf)};
    my @s=();
    my @e=();
    get_boundries(\@s,\@e,$Vec);
    my $gid=0;
    Output(1,$s[0]-1,$chr[$i],\$gid);
    for(my $j=0;$j<@s-1;$j++){
	Output($e[$j]+1,$s[$j+1]-1,$chr[$i],\$gid);
    }
    my $c="chr".$chr[$i];
    print STDERR get_seq_length($c,$gd),"\n";
    Output($e[@e-1]+1,get_seq_length($c,$gd),$chr[$i],\$gid);
    print STDERR "OK!\n";
}
system("rm -rf tmp");
exit;


########## sub routines ###########

sub get_seq_length{ 
    my ($chr, $dir) = @_;
    my $seq="";
    unless($dir=~/\/$/){
	$dir.="/";
    }
    open(SEQ,"$dir$chr.fa")||die(print "Open chrN.fa file error!\n");  #Open chrN.fa
    while(<SEQ>){ 
	chomp; 
	if (/^\>(\w+)/) {
	    $seq = ""; 
	}
	else{
	    $seq .= $_; 
	}
    }
    close(SEQ);  
    return length($seq);
}

sub Output{
    my ($s,$e,$chr,$gad)=@_;
    ${$gad}++;
    print "chr$chr\thg19_ref_gene\tintergenic\t$s\t$e\t0\t+\t.\tgene_id \"chr$chr.${$gad}\"; transcript_id \"chr$chr.${$gad}\";\n";
}

sub get_boundries{
    my ($sad,$ead,$vec)=@_;
    foreach my $g (@{$vec}){
	my $s=${${${$g}[0]}[0]}{start};
	my $e=${${${$g}[0]}[@{${$g}[0]}-1]}{end};
	foreach my $t (@{$g}){
	    if($s<${${$t}[0]}{start}){
		$s=${${$t}[0]}{start};
	    }
	    if($e>${${$t}[@{$t}-1]}{end}){
		$e=${${$t}[@{$t}-1]}{end};
	    }
	}
	push @{$sad},$s;
	push @{$ead},$e;
    }
}
sub Read_GTF{
    my ($file) = @_;
    my @gvec;
    my @re;
    my $begin = 1;
    my $pre_tid;
    my $pre_gid;
    my $i = 0;
    my $j = 0;
    my $k = 0;
    my @temp;
    open(IN,$file)||die "open $file error:$!\n";
    while(<IN>){
	chomp;
	if(/^[0-9a-z]+/){

	    if($begin == 1){
		@temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[\w\.]+\"); (transcript_id \"[\w\._]+\")/;
		$pre_tid = $2;
		$pre_gid = $1;

		$gvec[$i] = init_gene();
		${$gvec[$i]}[$j] = init_tran();

		${${$gvec[$i]}[$j]}[$k] =
			  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		$begin = 0;
	    }
	    else{
		@temp = split("\t",$_);
		$temp[8] =~ /(gene_id \"[\w\.]+\"); (transcript_id \"[\w\._]+\")/;

		if($1 eq $pre_gid){
		    if($2 eq $pre_tid){
			$k++;
			${${$gvec[$i]}[$j]}[$k] =
				  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		    }
		    else{
			$j++;
			${$gvec[$i]}[$j] = init_tran();
			$k = 0;
			${${$gvec[$i]}[$j]}[$k] = 
				  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
			$pre_tid = $2;
		    }
		}
		else{
		    $i++;
		    $gvec[$i] = init_gene();
		    $j = 0;
		    ${$gvec[$i]}[$j] = init_tran();
		    $k = 0;
		    ${${$gvec[$i]}[$j]}[$k] = 
			      new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		    $pre_gid =$1;
		    $pre_tid = $2;
		}
	    }
	}
    }
    close IN;
    $re[0] = \@gvec;
    $re[1] = scalar(@gvec);
    return \@re;
}

sub new_sf{             
    my ($chr,$src,$type,$start,$end,$score,$sym,$phase,$gid,$tid)=@_;
    my %hash = (
	"chr" => $chr,
	"source" => $src,
	"type" => $type,
	"score" => $score,
	"gid" => $gid,
	"tid" => $tid,
	"start" => $start,
	"end" => $end,
	"sym" => $sym,
	"phase" => $phase,
    );
    return \%hash;
}

sub init_gene{
    my @a;
    return \@a;
}
sub init_tran{
    my @a;
    return \@a;
}
