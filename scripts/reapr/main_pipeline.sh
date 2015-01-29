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

rhody=true
MP_lib=false
PE_lib=false
cov=100

reapr=true
getdistr=true


if [ "$rhody" = true ] ; 
then
        rm -r 'data-'"$rhody"'-ref'
        mkdir 'data-'"$rhody"'-ref'
        out='data-'"$rhody"'-ref'
        BASEDIR='/proj/b2013169/private/data/getdistr/reapr/'"$rhody"'-ref/'
        ref='/proj/b2013169/private/data/bacgenomes/references/modified_references/rhody.reference.fasta_modified_lines.fa'
        reads1='/proj/b2013169/private/data/bacgenomes/data/rhody_data/jump_reads.Solexa-62923.A.PE_orient.fasta'
        reads2='/proj/b2013169/private/data/bacgenomes/data/rhody_data/jump_reads.Solexa-62923.B.PE_orient.fasta'
        # input="$BASEDIR"'input/'
        # reapr_out="$BASEDIR"'reapr/' 
        # getdistr_out="$BASEDIR"'getdistr/' 

        window_size=2570
        declare -a gaps=("0" "500" "1000" "1500" "2000" "2500" "3000" "3500" "4000")
        declare -a errors=("-2000" "-1500" "-1000" "-500" "0" "500" "1000" "1500" "2000" )

elif [ "$MP_lib" = true ]

        mean=3000
        stddev=500
        min=1500
        max=4500
        window_size=2801
        declare -a gaps=("0" "500" "1000" "1500" "2000")
        declare -a errors=("-2000" "-1500" "-1000" "-500" "0" "500" "1000" "1500" "2000" )

        rm -r 'data-'"$mean-$stddev-$cov"
        mkdir 'data-'"$mean-$stddev-$cov"
        out='data-'"$mean-$stddev-$cov"
        BASEDIR='/proj/b2013169/private/data/getdistr/reapr/'"$mean-$stddev-$cov/"
        # input="$BASEDIR"'input/' 
        # reapr_out="$BASEDIR"'reapr/' 
        # getdistr_out="$BASEDIR"'getdistr/' 


else
        mean=500
        stddev=50
        min=350
        max=650
        window_size=300
        declare -a gaps=("0" "50" "100" "150" "200")
        declare -a errors=("-200" "-150" "-100" "-50" "0" "50" "100" "150" "200" )

        rm -r 'data-'"$mean-$stddev-$cov"
        mkdir 'data-'"$mean-$stddev-$cov"
        out='data-'"$mean-$stddev-$cov"
        BASEDIR='/proj/b2013169/private/data/getdistr/reapr/'"$mean-$stddev-$cov/"
        # input="$BASEDIR"'input/'
        # reapr_out="$BASEDIR"'reapr/' 
        # getdistr_out="$BASEDIR"'getdistr/' 
fi

if [ ! -d $BASEDIR ]; then
        mkdir -p $BASEDIR
fi


# if [ ! -d $reapr_out ]; then
#         mkdir $reapr_out
# fi
# if [ ! -d $getdistr_out ]; then
#         mkdir $getdistr_out
# fi

echo "Simulating"

for gap in "${gaps[@]}"
do
        if [ "$rhody" = true ] ; 
        then
                for error in "${errors[@]}"
                do
                        input="$BASEDIR"'input/gap_'"$gap"'_error_'"$error"
                        if [ ! -d $input ]; then
                                mkdir -p $input
                        fi
                        python /home/kris/git_repos/genomics_tools/scripts/rhodo_sim_individual_gaps.py \
                                "$ref"  "$reads1" "$reads2" "$input" "$gap" 
                done
        else

                input_unif="$BASEDIR"'input/normal_gap_'"$gap"
                input_normal="$BASEDIR"'input/uniform_gap_'"$gap"
                if [ ! -d $input_unif ]; then
                        mkdir -p $input_unif
                fi
                if [ ! -d $input_normal ]; then
                        mkdir -p $input_normal
                fi

                python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py \
                10000 $gap normal "$cov" 100 "$input_normal" -sort -scafs -errors "${errors[@]}"  -burnin 4000000 -nrgaps 100 --mean "$mean" --sd "$stddev" 
                
                python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py \
                10000 $gap uniform "$cov" 100 "$input_unif" -sort -scafs -errors "${errors[@]}"  -burnin 4000000 -nrgaps 100 --min_size "$min" --max_size "$max"

        fi
       #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 5000 0 300 30 10 50 "$input"'gap_'"$gap" -sort -scafs -errors -1500 -1250 -1000 -750 -500 -250 0 250 500 750 1000 1250 1500 -burnin 4000 -nrgaps 10
        #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 7000 $gap 300 30 100 50 "$input"'gap_'"$gap/$index" -sort -scafs -errors "$error" -burnin 40000 -nrgaps 10
        #python /home/kris/git_repos/genomics_tools/scripts/main_simulate_reapr.py 10000 $gap 2000 500 100 100 "$input"'gap_'"$gap/$index" -sort -scafs -errors "$error" -burnin 4000000 -nrgaps 100
        
done




# running reapr

if [ "$rhody" = true ] ; 
then
        for assembly in reference
        do

                for error in "${errors[@]}"
                do
                        reapr_input="$BASEDIR"'input/gap_'"$gap"'_error_'"$error"'/mapped.bam'
                        reapr_output="$BASEDIR"'reapr/'"$assembly"'/gap_'"$gap"'_error_'"$error"

                        if [ ! -d $reapr_out ]; then
                                mkdir -p $reapr_out
                        fi

                        echo "running reapr" "$reapr_out"

                        rm -r "$reapr_out"
                        reapr pipeline "$ref" "$reapr_input"  "$reapr_output"
                done
        done
else

        for distr in normal uniform
        do
                for gap in "${gaps[@]}"
                do
                        reapr_input="$BASEDIR"'input/'"$distr"'_gap_'"$gap"
                        reapr_output="$BASEDIR"'reapr/'"$distr"'/gap_'"$gap"
                        if [ ! -d $reapr_out ]; then
                                mkdir -p $reapr_out
                        fi                        
                        rm -r "$reapr_output" 

                        echo "running reapr" "$reapr_output"
                        reapr pipeline "$reapr_input"'/ctgs.fa' "$reapr_input"'/mapped.bam' "$reapr_output"
                
                done
        done

fi





# PARSE reaprs results

if [ "$rhody" = true ] ; 
then

        for assembly in reference
        do

                reapr_results="$out"'/reapr-'"$assembly"'.txt'


                if [ ! -f "$reapr_out$assembly"'/04.break.broken_assembly.fa' ]; then
                   echo "Reapr did not find any errors"
                   python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$input"'true_error_pos.gff' "$reapr_out$assembly"'/03.score.errors.gff' "$input"'ctgs.fa' >> "$reapr_results"
                   continue
                fi
                  python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$input"'true_error_pos.gff' "$reapr_out$assembly"'/03.score.errors.gff.gz' "$input"'ctgs.fa' >> "$reapr_results"
        done

else

        for distr in normal uniform
        do 
                
                #mkdir 'data'"$d"
                reapr_results="$out"'/reapr-'"$distr"'-.txt'
                touch "$reapr_results"
                for gap in "${gaps[@]}"
                do

                        if [ ! -f "$reapr_out$distr"'_gap_'"$gap/04.break.broken_assembly.fa" ]; then
                                echo "Reapr did not find any errors"
                                python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$input$distr"'_gap_'"$gap/true_error_pos.gff" "$reapr_out$distr"'_gap_'"$gap/03.score.errors.gff" "$input$distr"'_gap_'"$gap/ctgs.fa" >> "$reapr_results"
                                continue
                        fi 
                        python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py "$input$distr"'_gap_'"$gap/true_error_pos.gff" "$reapr_out$distr"'_gap_'"$gap/03.score.errors.gff.gz" "$input$distr"'_gap_'"$gap/ctgs.fa" >> "$reapr_results"
                        #python /home/kris/git_repos/GetDistr/scripts/reapr/parse_reapr_out.py  "$reapr_out"'gap_'"$gap/05.summary.stats.tsv" "reapr_results.txt"
                done
        done

fi



# Running our detector
echo "running KS detector" 

if [ "$rhody" = true ] ; 
then

        echo "running KS detector"
        for assembly in reference
        do
                rm -r "$getdistr_out$assembly"
                mkdir "$getdistr_out$assembly"
                echo "running KS" "$getdistr_out$assembly"
                python /home/kris/git_repos/GetDistr/scripts/reapr/find_misassemblies.py \
                python /home/kris/git_repos/GetDistr/scripts/reapr/find_misassemblies.py pipeline \
                "$input"'mapped.bam' "$ref" "$getdistr_out$assembly"  "$window_size"
        done

else

        for distr in normal uniform
        do 

                for gap in "${gaps[@]}"
                do
                        rm -r "$getdistr_out$distr"'_gap_'"$gap"
                        mkdir "$getdistr_out$distr"'_gap_'"$gap"
                        echo "running KS" "$getdistr_out$distr"'_gap_'"$gap/$error"
                        python /home/kris/git_repos/GetDistr/scripts/reapr/find_misassemblies.py \
                                "$input$distr"'_gap_'"$gap/mapped.bam"  \
                                "$input$distr"'_gap_'"$gap/ctgs.fa" \
                                "$window_size" \
                                "$getdistr_out$distr"'_gap_'"$gap"

                done
        done

fi



#Parse our results

if [ "$rhody" = true ] ; 
then

        for assembly in reference
        do
                getdistr_results="$out"'/getdistr-'"$assembly"'.txt'
                python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction_rhody.py \
                "$input"'true_error_pos.gff' \
                "$getdistr_out$assembly"'/estimated_misassm.gff' \
                "$input"'ctgs.fa' \
                >> "$getdistr_results"
        done
else

for distr in normal uniform
do 

        getdistr_results="$out"'/getdistr-'"$distr"'-.txt'
        touch "$getdistr_results"
        for gap in "${gaps[@]}"
        do
                python /home/kris/git_repos/GetDistr/scripts/reapr/parse_assembly_correction.py \
                "$input$distr"'_gap_'"$gap/true_error_pos.gff" \
                "$getdistr_out$distr"'_gap_'"$gap/estimated_misassm.gff" \
                "$input$distr"'_gap_'"$gap/ctgs.fa"\
                >> "$getdistr_results"
                
        done
done

fi


# PLOT RESULTS

if [ "$rhody" = true ] ; 
then
        for assembly in reference
        do
                python /home/kris/git_repos/GetDistr/scripts/reapr/plot_fp_tptry2.py \
                        "$reapr_results" "$getdistr_results" "$out"'/results_rhody_'"$assembly"'_'"$gap"
        done
else

        for distr in normal uniform
        do 

                getdistr_results="$out"'/getdistr-'"$distr"'-.txt'
                reapr_results="$out"'/reapr-'"$distr"'-.txt'

                for gap in "${gaps[@]}"
                do
                        
                        python /home/kris/git_repos/GetDistr/scripts/reapr/plot_fp_tptry2.py \
                                "$reapr_results" "$getdistr_results" "$out"'/results_'"$distr"'_'"$gap"

                done
        done

fi



# run QUAST
# for gap in 0 #250  500 750 1000 1250 1500
# do
#         for error in 1 # {1..13};
#         do
#                 if [ ! -f "$reapr_out"'gap_'"$gap/$error/04.break.broken_assembly.fa" ]; then
#                         echo "Reapr did not find any errors QUASTing original assembly"
#                         QUAST -s -R "$input"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'reapr/gap_'"$gap/$error"  "$input"'gap_'"$gap/$error"/ctgs.fa
#                         continue
#                 fi 
#                 QUAST -s --no-plots -R "$input"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'reapr/gap_'"$gap/$error"  "$reapr_out"'gap_'"$gap/$error"/04.break.broken_assembly.fa
#                 QUAST -s --no-plots -R "$input"'gap_'"$gap/$error"/genome.fa -o "$quast_out"'original/gap_'"$gap/$error"  "$input"'gap_'"$gap/$error"/ctgs.fa
#         done
# done

# parsing QUAST


# Latex table





