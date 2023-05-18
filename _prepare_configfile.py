#!/usr/bin/env python3
#

from os import listdir
import yaml
import pandas as pd


config = {
    "WDIR": "results",
}


def main():
    df = pd.read_csv("samples.txt", sep="\t", names=["r1", "r2", "sample_name"])

    samples = sorted(df[~df["sample_name"].str.contains("Neg")].sample_name.unique())

    config["SAMPLES"] = samples

    with open("config.yaml", 'w') as fd_config:
        yaml.dump(config, fd_config)


if __name__ == "__main__":
    main()
