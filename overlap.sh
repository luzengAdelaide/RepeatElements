#!/bin/bash

# analyze regulatory elements with repeats in gene regions
cd /scratch/luAnalysis/work_2013/regulator_gene_distribution

arr=('Blood cancer' 'Blood normal' 'Liver cancer' 'Lung normal' 'Bresat normal' 'Blood vessel normal')
arr1=('active_pro' 'weak_pro' 'strong_enhan' 'weak_enhan' 'repressed' 'insulator')
arr2=('5UTR' 'CDS_exon' 'CDS_intron' '3UTR' 'Intergenic')

for i in "${arr[@]}"
do

    for func in "${arr2[@]}"
    do

	for features in "${arr3[@]}"
	do
	    
	    echo $i $func $features
	    
	    intersectBed -a "$i"_"$func" -b Repeatmasker > repeat_"$i"_"func"
	    intersectBed -a "$func"_data -b repeat_"$i"_"func" > repeat_"func"_"$i"_"func"

	done
    done	
done