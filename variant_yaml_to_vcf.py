#!/usr/bin/env python

import yaml
import pyvcf
import argparse
from pathlib import Path


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def parse_yaml(yaml_path: Path) -> dict:
    """
    Read and parse the variant yaml into a dictionary
    """
    with open(yaml_path) as fh:
        yaml_variant = yaml.load(fh)
    return yaml_variant


def generate_vcf(variant_def: dict) -> vcf.Reader:
    """

    """





if __name__ == '__main__':

    parser = argparse.ArgumentParser("Convert PHE ncov yaml variant definitions to VCF format")
    parser.add_argument("-i", "--input", required=True, type=check_file,
                         help="Path to PHE genomics variant definition YAML file")

    # if output name isn't supplied then
    #parser.add_argument("-o", "--output", default=

    args = parser.parse_args()

    parse_yaml(args.input)


