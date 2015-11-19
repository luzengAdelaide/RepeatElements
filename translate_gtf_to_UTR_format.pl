#!/usr/bin/perl -w
use strict;


my ($gtf) = @ARGV;
my $usage = "$0 <gtf_file> 
Add 5_UTR and 3_UTR information into genes
";
die $usage if @ARGV<1;

my $utr_file = "UTR_format." .$gtf;
open(OUT, "> $utr_file") || die "can not open the $utr_file\n";

my $ad = read_gtf_to_tvec($gtf);
foreach my $t (@{$ad}){
    my $st=0;
    my $en=0;
    foreach my $i (@{$t}){
	if(${$i}{type} eq "start_codon"){
	    if(${$i}{sym} eq '+'){
		$st=${$i}{start};
	    }
	    else{
		$st=${$i}{end};
	    }
	}
	elsif(${$i}{type} eq "stop_codon"){
	    if(${$i}{sym} eq '+'){
		$en=${$i}{end};
	    }
	    else{
		$en=${$i}{start};
	    }
	}
    }
    if(${${$t}[0]}{sym} eq '+'){
	foreach my $i (@{$t}){
	    if(${$i}{type} eq "exon"){
		if(${$i}{end}<$st){
		    Output($i,${$i}{start},${$i}{end},"5_UTR");
		}
		elsif(${$i}{start}>$en){
		    Output($i,${$i}{start},${$i}{end},"3_UTR");
		}
		elsif(${$i}{start}<=$st && ${$i}{end}>=$en){
		    if(${$i}{start}<$st && ${$i}{end}>$en){
			Output($i,${$i}{start},$st-1,"5_UTR");
			Output($i,$en+1,${$i}{end},"3_UTR");
		    }
		    elsif(${$i}{start}<=$st && ${$i}{end}>$en){
			Output($i,$en+1,${$i}{end},"3_UTR");
		    }
		    elsif(${$i}{start}<$st && ${$i}{end}>=$en){
			Output($i,${$i}{start},$st-1,"5_UTR");
		    }
		}
		elsif(${$i}{start}<$st && ${$i}{end}>$st){
		    Output($i,${$i}{start},$st-1,"5_UTR");
		}
		
		elsif(${$i}{start}<$en && ${$i}{end}>$en){
		    Output($i,$en+1,${$i}{end},"3_UTR");
		}
	    }
	    else{
		Output($i,${$i}{start},${$i}{end},${$i}{type});
	    }
	}
    }
	
    else{
	foreach my $i (@{$t}){
	    if(${$i}{type} eq "exon"){
		if(${$i}{end}<$en){
		    Output($i,${$i}{start},${$i}{end},"3_UTR");
		}
		elsif(${$i}{start}>$st){
		    Output($i,${$i}{start},${$i}{end},"5_UTR");
		}
		elsif(${$i}{start}<=$en && ${$i}{end}>=$st){
		    if(${$i}{start}<$en && ${$i}{end}>$st){
			Output($i,${$i}{start},$en-1,"3_UTR");
			Output($i,$st+1,${$i}{end},"5_UTR");
		    }
		    elsif(${$i}{start}<=$en && ${$i}{end}>$st){
			Output($i,$st+1,${$i}{end},"5_UTR");
		    }
		    elsif(${$i}{start}<$en && ${$i}{end}>=$st){
			Output($i,${$i}{start},$en-1,"3_UTR");

		    }
		}
		elsif(${$i}{start}<$en && ${$i}{end}>$en){
		    Output($i,${$i}{start},$en-1,"3_UTR");
		}
		
		elsif(${$i}{start}<$st && ${$i}{end}>$st){
		    Output($i,$st+1,${$i}{end},"5_UTR");
		}
	    }
	    else{
		Output($i,${$i}{start},${$i}{end},${$i}{type});
	    }
	}
    }
}
############### sub_routines #################
sub Output{
    my $i=$_[0];
    print OUT "${$i}{chr}\t${$i}{source}\t$_[3]\t$_[1]\t$_[2]\t${$i}{score}\t${$i}{sym}\t${$i}{phase}\t${$i}{gid}; ${$i}{tid};\n";
}
sub read_gtf_to_tvec{
    my $i = 0;
    my @tvec;
    open(A,"$_[0]");
    my $flag = 0; #flag = 0 init_state
    while(<A>){
	chomp;
	if($flag ==1){
	    if(/^[0-9a-z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[\w\_]+\"); (transcript_id \"[\w\_]+\")/;
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
	    if(/^[0-9a-z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[\w\.]+\"); (transcript_id \"[\w\.]+\")/;
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
