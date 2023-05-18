#!/usr/bin/env python3
#

import pandas as pd
from subprocess import check_call
from os import listdir
import json


def main():
    folder = "/mnt/data/bigdata/ideARG_2020/delivery_20230216/raw_reads"

    df = pd.read_csv("samples.txt", sep="\t", names=["r1", "r2", "sample_name"])
    for i in df.itertuples():
        smp = i.sample_name
        if "Neg" in smp:
            continue
        
        r1 = i.r1
        r2 = i.r2
        assert r1 in listdir(folder), f"R1 {r1}\n{listdir(folder)}"
        assert r2 in listdir(folder), f"R2 {r2}\n{listdir(folder)}"
        folder_sample = f"results/{smp}"
        check_call(f"mkdir -p {folder_sample}", shell=True)
        check_call(f"ln -s -T {folder}/{r1} {folder_sample}/R1.fastq.gz", shell=True)
        check_call(f"ln -s -T {folder}/{r2} {folder_sample}/R2.fastq.gz", shell=True)


if __name__ == "__main__":
    main()
