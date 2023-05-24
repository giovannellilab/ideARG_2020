#!/bin/bash

grep "^NAME" Resfams-full.hmm | awk '{print $2}' > hmms_key_names.txt

mkdir -p hmms_folder

while read name
do
    hmmfetch Resfams-full.hmm $name > hmms_folder/$name.hmm
done < hmms_key_names.txt
