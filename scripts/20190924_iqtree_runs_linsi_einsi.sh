#!/bin/bash

#runs iqtree with a given model N_runs times 
#
#Current possible folders to run it on are:
#
#GPM2_GPM3_OG6006  :   best model: LG+F+I+G4
#
#Usage:
#/home/heineike/github/y1000plus_tools/scripts/20190924_iqtree_runs.sh <base_dir> <N_runs> <outgroup_id> <model>
# 
# -bnni flag does NNI search during bootstrap - takes longer but more thorough for "severe model violations"
# -wbtl prints out the bootstrap trees


base_dir=/home/heineike/genomes/y1000plus/proteins_og/

aln_dir=${base_dir}$1

N_runs=$2

#Activate python environment for iqtree
. activate /home/lab/envs/seqanalysis

mkdir ${aln_dir}/tree

for jj in $(seq 1 $N_runs)
do
    echo finding models run $jj
    aln_base_name_einsi=${aln_dir}/${1}_aln_einsi_trimmed.fasta
    aln_base_name_linsi=${aln_dir}/${1}_aln_linsi_trimmed.fasta
    iqtree -s ${aln_base_name_einsi} -nt AUTO -o $3 -m $4 -bb 1000 -alrt 1000 -bnni -wbtl
    iqtree -s ${aln_base_name_linsi} -nt AUTO -o $3 -m $4 -bb 1000 -alrt 1000 -bnni -wbtl
    mkdir ${aln_dir}/tree/run_einsi_${jj}
    mkdir ${aln_dir}/tree/run_linsi_${jj}
    mv ${aln_base_name_einsi}.* ${aln_dir}/tree/run_einsi_${jj}/
    mv ${aln_base_name_linsi}.* ${aln_dir}/tree/run_linsi_${jj}/
done
