#!/usr/bin/env python3
#

import yaml
from Bio.SeqIO.FastaIO import SimpleFastaParser
from os import listdir


resfinder_files = [
    "aminoglycoside.fsa",
    "beta-lactam.fsa",
    "colistin.fsa",
    "fosfomycin.fsa",
    "fusidicacid.fsa",
    "glycopeptide.fsa",
    "macrolide.fsa",
    "misc.fsa",
    "nitroimidazole.fsa",
    "oxazolidinone.fsa",
    "phenicol.fsa",
    "pseudomonicacid.fsa",
    "quinolone.fsa",
    "rifampicin.fsa",
    "sulphonamide.fsa",
    "tetracycline.fsa",
    "trimethoprim.fsa",
]



def main():
    folder = "databases/resfinder_db"

    with open("databases/all_resfinder_db.fasta", "wt") as fo:    
        for f in resfinder_files:
            category = f.split(".fsa")[0]
            with open(f"{folder}/{f}") as fd:
                for idx, seq in SimpleFastaParser(fd):
                    assert "|" not in idx
                    fo.write(f">{idx}|{category}\n{seq}\n")


if __name__ == "__main__":
    main()
