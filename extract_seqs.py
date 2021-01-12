#!/usr/bin/env python

import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")

def filter_metadata(metadata_fp, query):
    """
    Parse and filter metadata
    """
    df = pd.read_csv(str(metadata_fp), sep='\t')
    df = df.query(query)
    return df


def filter_seqs(seqs_fp, filtered_metadata):
    """
    Parse and filter the seqs file based on the filtered metadata
    """
    seq_names = set(filtered_metadata['strain'].values)
    filtered_seqs = []
    for record in SeqIO.parse(str(seqs_fp), 'fasta'):
        if record.id in seq_names:
            filtered_seqs.append(record)
    return filtered_seqs


def write_output(filtered_metadata, filtered_seqs, output_prefix):
    """
    Write filtered results
    """
    filtered_metadata.to_csv(output_prefix + "_metadata.tsv", sep='\t',
                             index=False)

    with open(output_prefix + "_seqs.fasta", 'w') as out_fh:
        SeqIO.write(filtered_seqs, out_fh, 'fasta-2line')


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Extract sequences from nextfasta based on "
                                     "pandas filter on nextmeta")
    parser.add_argument('-m', '--nextmeta', type=check_file, required=True,
                        help="Path to nextmeta file")
    parser.add_argument('-f', '--nextfasta', type=check_file, required=True,
                        help="Path to nextfasta file")
    parser.add_argument("-q", "--query", type=str, required=True,
                         help="Pandas query e.g., \"pangolin_lineage=='B.1.1.28'\"")
    parser.add_argument("-o", "--output_prefix", default="filtered",
                        help="Prefix for output filtered metadata and seqs")

    args = parser.parse_args()

    filtered_metadata = filter_metadata(args.nextmeta, args.query)

    filtered_seqs = filter_seqs(args.nextfasta, filtered_metadata)

    write_output(filtered_metadata, filtered_seqs, args.output_prefix)
