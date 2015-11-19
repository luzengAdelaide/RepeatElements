#!/bin/bash

cd /scratch/luAnalysis/work_2013/regulator_gene_distribution

# calculate how many base pairs that overlapped between TEs, REs and gene regions


for i in *_25bp
do 

    perl count_bp.pl $i > $i_bp

done