#!/usr/bin/env python

import argparse
import subprocess
from pathlib import Path
import time
import pandas as pd
import pysam
import shutil


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def update_pangolin():
    """
    Ensure pangolin is updated to the latest release
    """
    subprocess.check_output(["pangolin", "--update"])


def run_pangolin(input_genomes, threads):
    """
    Execute pangolin and collect assignments
    """
    output_dir = Path(f"pangolin_tmp_{time.time()}")
    subprocess.check_output(f"pangolin {input_genomes} -t {threads} "
                            f"-o {str(output_dir)}".split(),
                            stderr=subprocess.DEVNULL)

    output_path = output_dir / "lineage_report.csv"
    if not output_path.exists():
        raise FileNotFoundError(f"{str(output_path)} not created, check "
                                 "pangolin install")

    pangolin_df = pd.read_csv(str(output_path), sep=',')

    # tidy up the dataframe
    pangolin_df = pangolin_df.rename(columns={'taxon': 'SPECIMEN_ID',
                                              'lineage': 'LINEAGE',
                                              'status': 'STATUS',
                                              'note': 'NOTE',
                                              'pangolin_version': 'PANGOLIN_VERSION',
                                              'pangoLEARN_version': 'PANGOLEARN_VERSION',
                                              'pango_version': 'PANGO_VERSION'})

    # select only reported columns
    pangolin_df = pangolin_df[['SPECIMEN_ID', 'LINEAGE', 'PANGOLIN_VERSION',
                               'PANGOLEARN_VERSION', 'PANGO_VERSION', 'STATUS',
                               'NOTE']]

    pangolin_df['SPECIMEN_ID'] = pangolin_df['SPECIMEN_ID'].astype(str)
    # remove temp output
    shutil.rmtree(output_dir)

    return pangolin_df


def genome_completeness(input_genomes, reference_genome, threads):
    """
    Just calculate genome completeness by %N from minimap2 alignment
    """
    output_dir = Path(f"minimap2_tmp_{time.time()}")
    output_dir.mkdir()
    subprocess.check_output(f"minimap2 -a -x asm5 --secondary=no "
                            f"-t {threads} "
                            f"{reference_genome} {input_genomes} "
                            f"-o {output_dir}/align.sam".split(),
                            stderr=subprocess.DEVNULL)

    output_path = output_dir / "align.sam"
    if not output_path.exists():
        raise FileNotFoundError(f"{str(output_path)} not created, check "
                                 "minimap2 install")

    completeness = {'SPECIMEN_ID': [],
                    'GENOME_COMPLETENESS': []}
    samfile = pysam.AlignmentFile(output_path, 'r')
    for genome in samfile.fetch():
        completeness['SPECIMEN_ID'].append(str(genome.query_name))
        completeness['GENOME_COMPLETENESS'].append(len(genome.get_aligned_pairs(matches_only=True)) / samfile.get_reference_length("MN908947.3") * 100)
    completeness = pd.DataFrame(completeness)

    # remove temp output
    shutil.rmtree(output_dir)
    return completeness


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Generate completeness and pangolin '
                                     'report in PHO format')
    parser.add_argument("-i", "--input_genomes", type=check_file, required=True,
                        help="Concatenated fasta containing consensus genomes")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Number of threads to use")
    parser.add_argument("-r", "--reference", type=check_file, required=True,
                        help="Path to MN908947.3 reference genome in fasta format")
    parser.add_argument("-o", "--output", default="pho_report.csv",
                        help="Output CSV containing reporting information")
    args = parser.parse_args()

    update_pangolin()
    pangolin_df = run_pangolin(args.input_genomes, args.threads)

    genome_df = genome_completeness(args.input_genomes, args.reference,
                                    args.threads)

    output = pd.merge(pangolin_df, genome_df, on='SPECIMEN_ID',
                      how='outer',
                      validate='one_to_one')

    output.to_csv(args.output, index=False, sep='\t')
