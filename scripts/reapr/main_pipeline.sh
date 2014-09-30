#!/bin/bash -l
#SBATCH -A b2013169
#SBATCH -p node -n 1
#SBATCH -t 1-5:00:00
#SBATCH -J REAPR

for arg in "$@"; do
    args="$args $arg"
done


#module add bioinfo-tools
#module add bwa

for gap in 0 250 500 750 1000 1250 1500 
do
	python /home/kris/git_repos/genomics_tools/scripts/main_simulate.py 10000 $gap 2000 500 100 100 /tmp/reapr_in -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000000 -nrgaps 100
   	echo output_reapr_$gap
   	reapr pipeline /tmp/reapr_in/ctgs.fa /tmp/reapr_in/mapped.bam "$args_$gap"
   	python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$args_$gap/05.summary.report.txt" "reapr_results.txt" 
done

