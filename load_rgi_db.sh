#!/bin/bash
#


rgi clean

# temp_folder="./card_db"
temp_folder=$1

mkdir -p $temp_folder

cd $temp_folder

wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir -p wildcard
tar -xjf wildcard_data.tar.bz2 -C wildcard
gunzip wildcard/*.gz


rgi card_annotation -i card.json > card_annotation.log 2>&1

temp_var=$(ls *_all.fasta | grep "_v" | head -n 1)

IFS="_" read b1 b2 version remaining <<< $temp_var
version_rgi_db=$(echo $version | cut -c2-)

rgi wildcard_annotation \
    -i wildcard \
    --card_json card.json \
    -v $version_rgi_db > wildcard_annotation.log 2>&1


rgi load \
  --card_json card.json \
  --debug \
  --card_annotation card_database_v${version_rgi_db}.fasta \
  --card_annotation_all_models card_database_v${version_rgi_db}_all.fasta \
  --wildcard_annotation wildcard_database_v${version_rgi_db}.fasta \
  --wildcard_annotation_all_models wildcard_database_v${version_rgi_db}_all.fasta \
  --wildcard_index wildcard/index-for-model-sequences.txt \
  --wildcard_version ${version_rgi_db} \
  --amr_kmers wildcard/all_amr_61mers.txt \
  --kmer_database wildcard/61_kmer_db.json \
  --kmer_size 61

mv localDB ../


