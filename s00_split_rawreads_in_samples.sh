#!/bin/bash

folder_rawreads="delivery_20230216/raw_reads"

cd $folder_rawreads

for file in *_R1*.gz; do
	prefix="${file%_R1*}"
	suffix="${file#*_R1}"	

	file_out_1=$(basename "${prefix}_R1${suffix}")
	file_out_2=$(basename "${prefix}_R2${suffix}")

    IFS="-" read first_part sample_name sample_num remaining <<< $file_out_1

    printf "$file_out_1\t$file_out_2\t${sample_name}_${sample_num}\n"

done > samples.txt
