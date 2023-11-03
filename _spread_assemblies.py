#!/usr/bin/env python3
#

import yaml
from subprocess import check_call


def main():
    folder = "from_kbase/assemblies/metaspades"
    out_folder = "results_metaspades"    

    with open("config.yaml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    samples = config["SAMPLES"]

    for s in samples:
        metaspades = f"{out_folder}/{s}/metaspades"
        check_call(f"mkdir -p {metaspades}", shell=True)
        check_call(f"cp {folder}/{s}.metaspades.assembly.fa {metaspades}/contigs.fasta", shell=True)


if __name__ == "__main__":
    main()
