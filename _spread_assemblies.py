#!/usr/bin/env python3
#

import yaml
from subprocess import check_call


def main():
    folder = "from_kbase/assemblies"
    out_folder = "results"
    

    with open("config.yaml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    samples = config["SAMPLES"]

    for s in samples:
        megahit = f"{out_folder}/{s}/megahit"
        check_call(f"mkdir -p {megahit}", shell=True)
        check_call(f"cp {folder}/{s}.megahit.assembly.FASTA/{s}.megahit.assembly.fa {megahit}/contigs.fasta", shell=True)


if __name__ == "__main__":
    main()
