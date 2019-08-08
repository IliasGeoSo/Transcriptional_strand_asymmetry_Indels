#!/bin/sh

echo "############################### This is the start of the qsub.sh ###############################"
echo "This is the array index"
echo $LSB_JOBINDEX

/nfs/users/nfs_i/igs/miniconda2/bin/python /nfs/compgen-04/team218/ilias/All_chr_hg19_TSA_non_overlap/scripts_indels_paper/indels_asymmetry/simulated_indels/control_gen.py $LSB_JOBINDEX 


