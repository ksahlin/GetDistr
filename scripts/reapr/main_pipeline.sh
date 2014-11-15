#!/bin/bash -l
#SBATCH -A b2013169
#SBATCH -p node -n 1
#SBATCH -t 3-5:00:00
#SBATCH -J REAPR

# for arg in "$@"; do
#     args="$args $arg"
# done

reapr_in='/proj/b2013169/private/data/getdistr/reapr/input/'
reapr_out='/proj/b2013169/private/data/getdistr/reapr/output/'

module add bioinfo-tools
module add bwa

#       python /home/kris/git_repos/genomics_tools/scripts/main_simulate.py 10000 $gap 2000 500 100 100 "$reapr_in" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000000 -nrgaps 100

echo "Simulating"

for gap in 0 250 500 750 1000 1250 1500
do
        python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 5000 0 300 30 10 50 "$reapr_in"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000 -nrgaps 10
        #python /home/kris/git_repos/genomics_tools/scripts/main_simulate.py 10000 $gap 2000 500 100 100 "$reapr_in"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000000 -nrgaps 100
done

for gap in 0 250 500 750 1000 1250 1500
do
        #python /home/kris/git_repos/genomics_tools/scripts/main_simulate.py 10000 $gap 2000 500 100 100 "$reapr_in"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000000 -nrgaps 100
        echo "running reapr" "$reapr_out"'gap_'"$gap"
        rm -r "$reapr_out"'gap_'"$gap"
        reapr pipeline "$reapr_in"'gap_'"$gap"/ctgs.fa "$reapr_in"'gap_'"$gap"/mapped.bam "$reapr_out"'gap_'"$gap"
        #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
done


rm reapr_results.txt
touch reapr_results.txt

for gap in 0 250  500 750 1000 1250 1500
do
        python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$reapr_in"'gap_'"$gap/true_error_pos.gff" "$reapr_out"'gap_'"$gap/03.score.errors.gff.gz" "/tmp/reaper_out" >> "reapr_results.txt"
        #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
done