#!/usr/bin/perl -w
use strict;

my ($gtf) = @ARGV;
my $usage = "usage:   $0 [option] <gtf_annotate> >new_gtf
Merge transcripts to genes----change gene_id in gtf files.
";
die $usage if @ARGV<1;

my @chr = 1..22;
$chr[22] = 'X';
$chr[23] = 'Y';
my $g_num = 0;
my $t_num = 0;

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
    print STDERR "chr$chr[$i] processing:\n";
    my $ann_gtf = "tmp/annotation/".$chr[$i];
    
    print STDERR "\tRead chr$chr[$i] ......";
    my $ann_ad = read_gtf_to_tvec($ann_gtf);
    print STDERR "OK!\n";
    my @ann_keys_s;
    my @ann_keys_e;
    my @pre_keys_s;
    my @pre_keys_e;
    print STDERR "\tGet keys for chr$chr[$i] ......";
    get_keys(\@ann_keys_s,\@ann_keys_e,$ann_ad);
    print STDERR "OK!\n";
    print STDERR "\tSort for chr$chr[$i] annotation......";
    MergeSort($ann_ad,\@ann_keys_s,\@ann_keys_e);
    print STDERR "OK!\n";

    print STDERR "\tCheck and print for chr$chr[$i] ......";
    Merge_And_Print($ann_ad,\@ann_keys_s,\@ann_keys_e,\$g_num,\$t_num);
    print STDERR "OK!\n";
}
print STDERR "Total gene number : $g_num\n";
print STDERR "Total transcript number : $t_num\n";
 
system("rm -rf tmp");
exit;


############### sub_routines #################

sub Merge_And_Print{
    my ($ann_ad,$ann_ks,$ann_ke,$gad,$tad)=@_;
    my $gid = 1;
    my $tid = 1;
    my $e=0;
    for(my $i=0;$i<@{$ann_ad};$i++){
	if($i==0){
	    $e=${$ann_ke}[0];
	    OutputTran(${$ann_ad}[$i],$gid,$tid);
	}
	else{
	    if($e<${$ann_ke}[$i-1]){
		$e=${$ann_ke}[$i-1];
	    }
	    if(${$ann_ks}[$i]<=$e){
		$tid++;
		OutputTran(${$ann_ad}[$i],$gid,$tid);
	    }
	    else{
		$gid++;
		$tid=1;
		OutputTran(${$ann_ad}[$i],$gid,$tid);
	    }
	}
    }
    ${$tad}+=@{$ann_ad};
    if(@{$ann_ad}>0){
	${$gad}+=$gid;
    }
}



sub OutputTran{              #\@tvec,$name ,output file
    my ($ad,$g,$t) = @_;
    my $d1;
    for($d1=0;$d1<@{$ad};$d1++){
	print  "${${$ad}[$d1]}{chr}\t";
	print  "${${$ad}[$d1]}{source}\t";
	print  "${${$ad}[$d1]}{type}\t";
	print  "${${$ad}[$d1]}{start}\t";
	print  "${${$ad}[$d1]}{end}\t";
	print  "${${$ad}[$d1]}{score}\t";
	print  "${${$ad}[$d1]}{sym}\t";
	print  "${${$ad}[$d1]}{phase}\t";
	print  "gene_id \"${${$ad}[$d1]}{chr}.$g\"; \t";
	print  "${${$ad}[$d1]}{tid};";
	print  "${${$ad}[$d1]}{flag}\n";    
    }
#    print "\n";
}


sub get_keys{
    my ($keys_s,$keys_e,$ad) = @_;
    for(my $i=0;$i<@{$ad};$i++){
	${$keys_s}[$i] = ${${${$ad}[$i]}[0]}{start};
	${$keys_e}[$i] = ${${${$ad}[$i]}[0]}{end};
	for(my $j=0;$j<@{${$ad}[$i]};$j++){
	    if(${$keys_s}[$i] > ${${${$ad}[$i]}[$j]}{start}){
		${$keys_s}[$i] = ${${${$ad}[$i]}[$j]}{start};
	    }
	    if(${$keys_e}[$i] < ${${${$ad}[$i]}[$j]}{end}){
		${$keys_e}[$i] = ${${${$ad}[$i]}[$j]}{end};
	    }
	}
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
	    if(/^[0-9a-z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[\w\.]+\"); (transcript_id \"[\w\.]+\")/;
		if($2 eq ${$tvec[$i]}[0]{tid}){
		    put_tran($tvec[$i],new_sf($1,$2,\@temp));
		}
		else{
		    $i++;
		    $tvec[$i] = init_tran(); 
		    put_tran($tvec[$i],new_sf($1,$2,\@temp));
		}
	    }
	}
	
	else{
	    $tvec[$i] = init_tran(); 
	    if(/^[0-9a-z]+/){
		my @temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[\w\.]+\"); (transcript_id \"[\w\.]+\")/;
		put_tran($tvec[$i],new_sf($1,$2,\@temp));
	    }
	    $flag = 1;
	}
    }  
    close A;
    return \@tvec;
}

sub new_sf{
    my ($gid,$tid,$ad)=@_;
    my $i=0;
    my $flag="\t";
    while(defined(${$ad}[9+$i])){
	$flag=$flag.${$ad}[9+$i]."\t";
	$i++;
    }
    my ($chr,$src,$type,$start,$end,$score,$sym,$phase)=@{$ad};
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
		"flag" => $flag,
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
