#!/bin/bash

#runs iqtree with a given model N_runs times 
#
#Current possible folders to run it on are:
#
#GPM2_GPM3_OG6006  :   best model: LG+F+I+G4
#
#Usage:
#/home/heineike/github/y1000plus_tools/scripts/20190924_iqtree_runs.sh <base_dir> <N_runs> <outgroup_id> <model>
# $1 is the base directory
# $2 is the number of runs
# $3 is the outgroup id
# $4 is the model we want to use


base_dir=
#/home/heineike/genomes/y1000plus/proteins_og/

aln_dir=${base_dir}$1

N_runs=$2

#Activate python environment for iqtree if necessary
#. activate /home/lab/envs/seqanalysis

mkdir ${aln_dir}/tree

#-bb is bootstrap with 1000 bootstrap iterations
#-bnni is the slower, but more thorough bootstrap search
#-alrt finds sh-alrt values, 1000 iterations


for jj in $(seq 1 $N_runs)
do
    echo finding models run $jj
    aln_base_name=${aln_dir}/${1}_aln_trimmed.fasta
    iqtree -s ${aln_base_name} -nt AUTO -o $3 -m $4 -bb 1000 -bnni -alrt 1000 
    mkdir ${aln_dir}/tree/run${jj}
    mv ${aln_base_name}.* ${aln_dir}/tree/run${jj}/
done
