#!/usr/bin/env python
from pathlib import Path
import pandas as pd
import upsetplot
from upset_plotly import plot
import argparse
import os


def parse_type_variant(input_file: Path) -> upsetplot.UpSet:
    """
    Parse the output csv from type_variant
    """
    variants = pd.read_csv(input_file, sep=',')

    variant_sets = {}
    genotype_columns = []
    for column in variants.columns:
        # handle SNP and amino acid changes
        if column.startswith('aa:') or column.startswith('snp:'):
            variant = column
            # last character of variant
            alt = variant[-1]
            strains_with_alt = set(variants.loc[variants[variant]==alt,
                                                'query'].values)
            variant_sets[variant] = strains_with_alt
        # handle deletions slightly differently
        elif column.startswith('del:'):
            variant = column
            strains_with_del =  set(variants.loc[variants[variant]=='del',
                                                 'query'].values)
            variant_sets[variant] = strains_with_del

            # x is returned if something other than ref or deletion is found
            strains_with_other = set(variants.loc[variants[variant]=='X',
                                                  'query'].values)
            if len(strains_with_other) > 0:
                variant_sets[variant + ":X"] = strains_with_other

        # skip other non-genotype columns
        else:
            pass

    # convert sets to df containing intersections
    variant_sets = upsetplot.from_contents(variant_sets)
    return variant_sets


def parse_ncov_watch(input_file: Path) -> upsetplot.UpSet:
    """
    Parse output from ncov_watch
    """
    variants = pd.read_csv(input_file, sep='\t')
    variant_sets = {}
    variant_sets = variants.groupby('mutation')['sample'].apply(list).to_dict()
    variant_sets = upsetplot.from_contents(variant_sets)
    return variant_sets


def variant_intersections(input_file: Path, input_type: str, output_path: Path):
    """
    Parse input and generate a plot using the upset_plotly library
    """

    if input_type == 'type_variants':
        variant_sets = parse_type_variant(input_file)
    elif input_type == 'ncov_watch':
        variant_sets = parse_ncov_watch(input_file)
    else:
        raise ValueError("input_type must be 'type_variants', "
                         f"'ncov_watch': {input_type}")

    fig = plot.upset_plotly(variant_sets, "Shared Variants Across Samples")

    variant_sets.to_csv(str(output_path) + ".tsv", sep='\t')

    if output_path.suffix == '.html':
        fig.write_html(str(output_path))
    elif output_path.suffix == '.png':
        fig.write_image(str(output_path))


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def check_output_suffix(path: str) -> Path:
    """
    Check specified output has .html or .png
    """
    path = Path(path)
    if path.suffix in ['.html', '.png']:
        return path
        #if os.access(path, os.W_OK):
        #    return path
        #else:
        #    raise argparse.ArgumentTypeError(f"{path} cannot be written")
    else:
        raise argparse.ArgumentTypeError(f"{path} must end in .html or .png")


if __name__ == '__main__':

    """Console script for variant_intersections"""
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', "--input", required=True, type=check_file,
                         help="Input type_variants or ncov_watch file")
    parser.add_argument('-t', '--type', required=True,
                         choices=['type_variants', 'ncov_watch'],
                         help="Program used to generate the input file")
    parser.add_argument('-o', '--output', default="output.html",
                        type=check_output_suffix,
                        help="Output file name (ending in .html for an "
                             "interactive plot or .png for static image)")

    args = parser.parse_args()
    variant_intersections(args.input, args.type,
                          args.output)

