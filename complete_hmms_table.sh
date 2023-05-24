#!/bin/bash
#

head -n1 results/cc4_1/hmms_presence/presence_results.tsv > complete_hmm_presence.tsv
awk -F '\t' 'FNR>1' ./results/**/hmms_presence/presence_results.tsv >> complete_hmm_presence.tsv
