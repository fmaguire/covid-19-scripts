#!/usr/bin/env python

import pandas as pd
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


def extract_protein_changes(row, gene):
    """
    Extract protein changes from nextclade output tsv
    """
    changes = []
    if not pd.isna(row['aaDeletions']):
        for deletion in row['aaDeletions'].split(','):
            if gene == 'all':
                changes.append(deletion)
            else:
                if deletion.startswith(f"{gene}:"):
                    changes.append(deletion)

    if not pd.isna(row['aaSubstitutions']):
        for mut in row['aaSubstitutions'].split(','):
            if gene == 'all':
                changes.append(mut)
            else:
                if mut.startswith(f"{gene}:"):
                    changes.append(mut)

    # sort based on gene and position
    changes = sorted(changes, key=lambda x: (x[0], int(x.split(':')[1][1:-1])))
    # convert back to string for saving as a flat text file
    changes = ",".join(changes)
    return changes


def add_mutations_to_metadata(metadata, nextclade, gene):
    """
    Add mutations/deletions to metadata file
    """
    nextclade[f"{gene} mutation/deletion sets"] = nextclade.apply(lambda x: extract_protein_changes(x, gene), axis=1)
    metadata[f"{gene} mutation/deletion sets"] = nextclade[f"{gene} mutation/deletion sets"]
    return metadata


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Extract nextclade mutations/deletions")
    parser.add_argument("-n", "--nextclade_output", required=True,
                         help="Path to nextclade output tsv")
    parser.add_argument("-m", "--metadata", required=True, type=check_file,
                         help="Path to nextmeta formatted metadata")
    parser.add_argument("-o", "--output", default="metadata_with_mutations.tsv",
                         help="Output path for metadata tsv with mutations")
    parser.add_argument("-g", "--gene", default="all", help="Gene for which to"
                         " extract mutations/deletions")

    args = parser.parse_args()

    nextclade = pd.read_csv(args.nextclade_output, sep='\t').set_index('seqName')

    metadata =  pd.read_csv(args.metadata, sep='\t').set_index('strain')

    metadata = add_mutations_to_metadata(metadata, nextclade, args.gene)

    metadata.reset_index().to_csv(args.output, sep='\t', index=False)
