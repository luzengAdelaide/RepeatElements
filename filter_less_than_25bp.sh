#!/bin/bash

cd /scratch/luAnalysis/work_2013/regulator_gene_distribution

# delete the overlapped regions less than 25bp

for i in repeat_*
do
    
    perl no_less_than_25bp.pl $i > $i_25bp

done

