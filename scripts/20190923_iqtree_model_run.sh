#!/bin/bash

#Runs model N_runs times for input alignment
#
#Current possible folders to run it on are: 
#
#GPM2_GPM3_OG6006
#
#Usage: 
#/home/heineike/github/y1000plus_tools/scripts/20190923_iqtree_model_run.sh <base_dir> <N_runs> <outgroup_id>

base_dir=/home/heineike/genomes/y1000plus/proteins_og/

aln_dir=${base_dir}$1

N_runs=$2

#Activate python environment for iqtree
. activate /home/lab/envs/seqanalysis

mkdir ${aln_dir}/model

for jj in $(seq 1 $N_runs)
do
    echo finding models run $jj
    aln_base_name=${aln_dir}/${1}_aln_trimmed.fasta
    iqtree -s ${aln_base_name} -nt AUTO -o $3 -m MF
    mkdir ${aln_dir}/model/run${jj}
    mv ${aln_base_name}.* ${aln_dir}/model/run${jj}/
done
