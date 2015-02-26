#!/bin/bash -l
#SBATCH -A b2013169
#SBATCH -p core -n 1
#SBATCH -t 02:00:00
#SBATCH -J staph-MP


#####
#Module dependencies
#none
#####

#########
#TO SPECIFY
INBASE='/proj/b2013169/private/data/genomes/staph/'
OUTBASE='/home/kris/Work/GetDistr/two_bias/staph/'

#######
/proj/b2013169/private/data/genomes/staph

#######
rm -r "$OUTBASE"
mkdir "$OUTBASE"


/usr/bin/time -v python /home/kris/git_repos/GetDistr/getdistr/assemblymodule/main.py pipeline \
       "$INBASE"aligned/ref/gage_mp.bam \
        "$INBASE"/ref/genome.fasta \
        "$OUTBASE" \
        --plots 1> "$OUTBASE"run.stdout 2> "$OUTBASE"run.stdout





