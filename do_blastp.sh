#!/bin/bash

read -p 'Number of threads (or cores) to use: ' THREADS
echo "Threads selected $THREADS"
SPLIT_FASTAS="./tmp/*"
RESULTS="$(ls ./tmp_results/*)"
DB="./Competition_data/cafa5protein/Train/train_sequences.fasta"
EVALUE="1e-10"
MAX_RESULTS="100"

for f in $SPLIT_FASTAS 
do 
    echo "$(basename -- $f)"
    #echo "$RESULTS"
    res="$(basename -- $f)"
    if echo "$RESULTS" | grep -qw "$res"
    then
    echo "Result $res Exists"
    else
    echo "Running $res"
    echo "query acc.ver,subject acc.ver,Percent identity,alignment length,mismatches,gap opens,q.start,q.end,s.start,s.end,evalue,bit score" > "./tmp_results/$res"
    var="$(./ncbi-blast-2.14/bin/blastp -db "$DB" -query "./tmp/$res" -evalue "$EVALUE" -max_target_seqs "$MAX_RESULTS" -outfmt 6 -num_threads "$THREADS")"
    echo "$var" | tr "\t" "," >> "./tmp_results/$res"
    fi
done 

exit 0