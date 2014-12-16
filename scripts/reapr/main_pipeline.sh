#!/bin/bash -l
#SBATCH -A b2013169
#SBATCH -p node -n 1
#SBATCH -t 3-5:00:00
#SBATCH -J REAPR

# for arg in "$@"; do
#     args="$args $arg"
# done


#####
#Module dependencies
module add bioinfo-tools
module add bwa
#####

#########
#TO SPECIFY
mean=500
stddev=50
min=350
max=650
cov=100
declare -a gaps=("0" "250" "500" "750" "1000" "1250" "1500")
declare -a errors=("-1500" "-1250" "-1000" "-750" "-500" "-250" "0" "250" "500" "750" "1000" "1250" "1500")
#######

#######
# Paths to files being generated
mkdir 'data-'"$mean-$stddev-$cov"
out='data-'"$mean-$stddev-$cov"
outdir='/proj/b2013169/private/data/getdistr/reapr/'"$mean-$stddev-$cov/"
reapr_in="$outdir"'input/' #'/proj/b2013169/private/data/getdistr/reapr/input/'
reapr_out="$outdir"'output/' # '/proj/b2013169/private/data/getdistr/reapr/output/'
#quast_out='/proj/b2013169/private/data/getdistr/reapr/quast/'
getdistr_out="$outdir"'getdistr/' # '/proj/b2013169/private/data/getdistr/reapr/getdistr/'

echo "Simulating"

for gap in "${gaps[@]}"
do
        for index in 1  # {1..13};
        do         
                let "error = -1500 + ($index -1) * 250"
                echo $error
                #for distr in normal uniform
                #do
                #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 5000 0 300 30 10 50 "$reapr_in"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000 -nrgaps 10
                #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 7000 $gap 300 30 100 50 "$reapr_in"'gap_'"$gap/$index" -sort -scafs -errors "$error" -burnin 40000 -nrgaps 10
                #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 10000 $gap 2000 500 100 100 "$reapr_in"'gap_'"$gap/$index" -sort -scafs -errors "$error" -burnin 4000000 -nrgaps 100
                python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 10000 $gap normal "$cov" 100 "$reapr_in"'normal_gap_'"$gap/$index" -sort -scafs -errors "${errors[@]}"  -burnin 10000000 -nrgaps 100 --mean "$mean" --sd "$stddev" 
                python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 10000 $gap uniform "$cov" 100 "$reapr_in"'uniform_gap_'"$gap/$index" -sort -scafs -errors "${errors[@]}"  -burnin 10000000 -nrgaps 100 --min_size "$min" --max_size "$max"
                #done
        done  
done




# running reapr
for distr in normal uniform
do
        for gap in "${gaps[@]}"
        do
                rm -r "$reapr_out$distr"'_gap_'"$gap" 
                mkdir "$reapr_out$distr"'_gap_'"$gap"   
                for error in  1 # {1..13};
                do 
                                echo "running reapr" "$reapr_out$distr"'_gap_'"$gap/$error"
                                #rm -r "$reapr_out$distr"'_gap_'"$gap/$error"

                        #python /home/kris/git_repos/genomics_tools/scripts/main_simulate.py 10000 $gap 2000 500 100 100 "$reapr_in"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000000 -nrgaps 100
                                reapr pipeline "$reapr_in$distr"'_gap_'"$gap/$error/ctgs.fa" "$reapr_in$distr"'_gap_'"$gap/$error/mapped.bam" "$reapr_out$distr"'_gap_'"$gap/$error"
                        #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
                done
        done
done



# PARSE reaprs results

for distr in normal uniform
do 
        
        #mkdir 'data'"$d"
        reapr_results="$out"'/reapr-'"$distr"'-.txt'
        touch "$reapr_results"
        for gap in "${gaps[@]}"
        do
                for error in 1 # {1..13};
                do       

                                if [ ! -f "$reapr_out$distr"'_gap_'"$gap/$error/04.break.broken_assembly.fa" ]; then
                                        echo "Reapr did not find any errors"
                                        python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$reapr_in$distr"'_gap_'"$gap/$error/true_error_pos.gff" "$reapr_out$distr"'_gap_'"$gap/$error/03.score.errors.gff" "$reapr_in$distr"'_gap_'"$gap/$error/ctgs.fa" >> "$reapr_results"
                                        continue
                                fi 
                                python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$reapr_in$distr"'_gap_'"$gap/$error/true_error_pos.gff" "$reapr_out$distr"'_gap_'"$gap/$error/03.score.errors.gff.gz" "$reapr_in$distr"'_gap_'"$gap/$error/ctgs.fa" >> "$reapr_results"
                                #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
                done
        done
done

# Running our detector
echo "running KS detector" 
for distr in normal uniform
do 

        for gap in "${gaps[@]}"
        do
                rm -r "$getdistr_out$distr"'_gap_'"$gap"
                mkdir "$getdistr_out$distr"'_gap_'"$gap"
                for error in 1 # {1..13};
                do 
                        echo "running KS" "$getdistr_out$distr"'_gap_'"$gap/$error"
                        python /home/kris/git_repos/GetDistr/scripts/reapr/FCD_test.py "$reapr_in$distr"'_gap_'"$gap/$error/mapped.bam"  "$reapr_in$distr"'_gap_'"$gap/$error/ctgs.fa"  0.05 "$mean" "$getdistr_out$distr"'_gap_'"$gap/$error"
                        #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
                done
        done
done

#Parse our results

for distr in normal uniform
do 
        # d=$(date '+%y-%m-%d-%H_%M_%S')
        # getdistr_results='getdistr-'"$distr"'-'"$d"'.txt'
        # touch "data/$getdistr_results"

        getdistr_results="$out"'/getdistr-'"$distr"'-.txt'
        touch "$getdistr_results"
        for gap in "${gaps[@]}"
        do
                for error in 1 # {1..13};
                do       

                        python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$reapr_in$distr"'_gap_'"$gap/$error/true_error_pos.gff" "$getdistr_out$distr"'_gap_'"$gap/$error/estimated_misassm.gff" "$reapr_in$distr"'_gap_'"$gap/$error/ctgs.fa" >> "$getdistr_results"
                done
        done
done

# PLOT RESULTS

for distr in normal uniform
do 

        getdistr_results="$out"'/getdistr-'"$distr"'-.txt'
        reapr_results="$out"'/reapr-'"$distr"'-.txt'
   
        python /home/kris/git_repos/GetDistr/scripts/reapr/plot_fp_tp.py "$getdistr_results" "$out"'/getdistr-'"$distr"
        python /home/kris/git_repos/GetDistr/scripts/reapr/plot_fp_tp.py "$reapr_results" "$out"'/reapr-'"$distr"

done


# run QUAST
# for gap in 0 #250  500 750 1000 1250 1500
# do
#         for error in 1 # {1..13};
#         do
#                 if [ ! -f "$reapr_out"'gap_'"$gap/$error/04.break.broken_assembly.fa" ]; then
#                         echo "Reapr did not find any errors QUASTing original assembly"
#                         QUAST -s -R "$reapr_in"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'reapr/gap_'"$gap/$error"  "$reapr_in"'gap_'"$gap/$error"/ctgs.fa
#                         continue
#                 fi 
#                 QUAST -s --no-plots -R "$reapr_in"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'reapr/gap_'"$gap/$error"  "$reapr_out"'gap_'"$gap/$error"/04.break.broken_assembly.fa
#                 QUAST -s --no-plots -R "$reapr_in"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'original/gap_'"$gap/$error"  "$reapr_in"'gap_'"$gap/$error"/ctgs.fa
#         done
# done

# parsing QUAST


# Latex table





