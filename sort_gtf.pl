#!/usr/bin/perl -w
use strict;

my ($gtf) = @ARGV;
my $usage = "$0 <gtf_file> > <output_gtf> 
";
die $usage if @ARGV<1;

my $ad = read_gtf_to_tvec($gtf);
foreach my $k (@{$ad}){
    SortGene($k);
    OutputTran($k);
}


exit;


############### sub_routines #################

sub SortGene{
    my ($ad)=@_;
    for(my $i=0;$i<@{$ad};$i++){
	for(my $j=$i+1;$j<@{$ad};$j++){
	    if($ad->[$i]->{start}>$ad->[$j]->{start}){
		($ad->[$i],$ad->[$j])=($ad->[$j],$ad->[$i]);
	    }
	    elsif($ad->[$i]->{start}==$ad->[$j]->{start}){
		if($ad->[$i]->{end}>$ad->[$j]->{end}){
		    ($ad->[$i],$ad->[$j])=($ad->[$j],$ad->[$i]);
		}
		elsif($ad->[$i]->{end}==$ad->[$j]->{end}){
		    if($ad->[$i]->{sym} eq '+' && ($ad->[$j]->{type} eq "start_codon" || $ad->[$i]->{type} eq "stop_codon")){
			($ad->[$i],$ad->[$j])=($ad->[$j],$ad->[$i]);
		    }
		    elsif($ad->[$i]->{sym} eq '-' && ($ad->[$j]->{type} eq "stop_codon" || $ad->[$i]->{type} eq "start_codon")){
			($ad->[$i],$ad->[$j])=($ad->[$j],$ad->[$i]);
		    }
		}
	    }
	}
    }
}

sub OutputTran{              #\@tvec,$name ,output file
    my ($ad) = @_;
    my $d1;
    for($d1=0;$d1<@{$ad};$d1++){
	print "${${$ad}[$d1]}{chr}\t";
	print "${${$ad}[$d1]}{source}\t";
	print "${${$ad}[$d1]}{type}\t";
	print "${${$ad}[$d1]}{start}\t";
	print "${${$ad}[$d1]}{end}\t";
	print "${${$ad}[$d1]}{score}\t";
	print "${${$ad}[$d1]}{sym}\t";
	print "${${$ad}[$d1]}{phase}\t";
	print "${${$ad}[$d1]}{gid}; ";
	print "${${$ad}[$d1]}{tid};\n";
	
    }
}

sub read_gtf_to_tvec{
    my $i = 0;
    my @tvec;
    open(A,"$_[0]");
    my $flag = 0; #flag = 0 init_state
    while(<A>){
	chomp;
	if($flag ==1){
	    if(/^[0-9a-zA-Z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[^\"]+\"); (transcript_id \"[^\"]+\")/;
		if($2 eq ${$tvec[$i]}[0]{tid}){
		    put_tran($tvec[$i],new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2));
		}
		else{
		    $i++;
		    $tvec[$i] = init_tran(); 
		    put_tran($tvec[$i],new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2));
		}
	    }
	}
    
	else{
	    $tvec[$i] = init_tran(); 
	    if(/^[0-9a-zA-Z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[^\"]+\"); (transcript_id \"[^\"]+\")/;
		put_tran($tvec[$i],new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2));
	    }
	    $flag = 1;
	}
    }  
    close A;
    return \@tvec;
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

sub init_tran{
    my @a;
    return \@a;
}
sub put_tran{       #tran array address,sf         
    my($tran,$sf)=@_;
    my $l = @{$tran};
    ${$tran}[$l] = $sf;    
}


sub MergeSort{
    my $vec = $_[0];
    my $num_ad = $_[1];
    my $rela_ad = $_[2];
    my $last=scalar(@{$num_ad})-1;
    &MergeSort_($vec,$num_ad,0,$last,$rela_ad);
}
 
sub MergeSort_ {
    my ($vec,$ad,$first,$last,$rela_ad)=@_;
    my $num_mid=int (($last+$first)/2);

    if($last-$first>0){
	&MergeSort_($vec,$ad,$first,$num_mid,$rela_ad);
	&MergeSort_($vec,$ad,$num_mid+1,$last,$rela_ad);
        &Merge($vec,$ad,$first,$num_mid,$last,$rela_ad);
    }
}


sub Merge {
    my ( $vec, $ad, $first, $num_mid, $last, $rela_ad ) =@_;
    my $i;
    my $f1=$first;
    my $f2=$num_mid+1;
    my @a;
    my @b;
    my @c;
    for($i=0;$i<=$last-$first;$i++){
	
	if($f1<=$num_mid && $f2<=$last){
	    if(${$ad}[$f1]<${$ad}[$f2]){
		$a[$i]=${$ad}[$f1];
		$b[$i]=${$vec}[$f1];
		$c[$i]=${$rela_ad}[$f1];
		$f1++;
	    }
	    elsif(${$ad}[$f1]>${$ad}[$f2]){
		$a[$i]=${$ad}[$f2];
		$b[$i]=${$vec}[$f2];
		$c[$i]=${$rela_ad}[$f2];
		$f2++;
	    }
	    else{
		if(${$rela_ad}[$f1]<=${$rela_ad}[$f2]){
		    $a[$i]=${$ad}[$f1];
		    $b[$i]=${$vec}[$f1];
		    $c[$i]=${$rela_ad}[$f1];
		    $f1++;
		}
		else{
		    $a[$i]=${$ad}[$f2];
		    $b[$i]=${$vec}[$f2];
		    $c[$i]=${$rela_ad}[$f2];
		    $f2++;
		}
	    }
	}
	
	elsif($f1<=$num_mid && $f2>$last){
	    $a[$i]=${$ad}[$f1];
	    $b[$i]=${$vec}[$f1];
	    $c[$i]=${$rela_ad}[$f1];
	    $f1++;
	}
	elsif($f1>$num_mid && $f2<=$last){
	    $a[$i]=${$ad}[$f2];
	    $b[$i]=${$vec}[$f2];
	    $c[$i]=${$rela_ad}[$f2];
            $f2++;
	}

    }
    for($i=0; $i<= $last-$first; $i++){
	${$ad}[$first+$i]=$a[$i];
	${$vec}[$first+$i]=$b[$i];
	${$rela_ad}[$first+$i]=$c[$i];
    }
}
