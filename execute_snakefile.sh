#!/bin/bash
#

# python3 _prepare_configfile.py

snakemake --use-conda --cores 30 $@ -s Snakefile.smk
