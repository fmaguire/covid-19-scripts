#!/usr/bin/env python
import pandas as pd
import gffutils
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


def summarise_ivar_variants(ivar_file, reference_gff_db):
    """
    parse ivar variants file, add product names from GFF file,
    then translate the mutations using summarise_ivar_mutations
    """
    var_df = pd.read_csv(ivar_file, sep='\t')

    var_df['isolate'] = ivar_file
    var_df = var_df.set_index('isolate')

    var_df['CDS_Product'] = var_df['GFF_FEATURE'].apply(lambda gff_feat: ref_db[gff_feat].attributes['Parent'][0] if not pd.isnull(gff_feat) else 'non-coding')
    var_df['CDS_Product'] = var_df['CDS_Product'].str.replace('gene-', '')
    var_df['Mutations'] = var_df.apply(describe_ivar_mutations, axis=1)
    return var_df


def describe_ivar_mutations(row):
    """
    Summarise mutations from ivar and place into the dataframe
    Specifically translating non-synonymous changes into amino acids
    """
    # stolen from https://github.com/cov-ert/type_variants/blob/master/type_variants.py#L17
    CDS_dict = {"orf1ab": (266, 21555),
                "s":       (21563, 25384),
                "orf3a":   (25393, 26220),
                "e":       (26245, 26472),
                "m":       (26523, 27191),
                "orf6":    (27202, 27387),
                "orf7a":   (27394, 27759),
                "orf8":    (27894, 28259),
                "n":       (28274, 29533),
                "orf10" :  (29558, 29674)}
    if row['ALT'].startswith('-'):
        # -2 to remove - and the anchor ref base
        return f"del:{row['POS']}:{len(row['ALT']) - 2}"

    elif row['CDS_Product'] == 'non-coding' or row['REF_AA'] == row["ALT_AA"]:
        return f"snp:{row['CDS_Product']}:{row['REF']}{row['POS']}{row['ALT']}"

    elif row['REF_AA'] != row["ALT_AA"]:
        aa_pos = int((row['POS'] - CDS_dict[row['CDS_Product'].lower()][0]) / 3 + 1)
        return f"aa:{row['CDS_Product']}:{row['REF_AA']}{aa_pos}{row['ALT_AA']}"


def parse_inputs(variant_file_list, ref_db):
	"""
	Parse set of input ivar variant files
	"""
	results = []
	# remove any inadvertent duplicate input files
	for ivar_variant_file in set(variant_file_list):
		ivar_variants = summarise_ivar_variants(ivar_variant_file, ref_db)
		results.append(ivar_variants)
	results = pd.concat(results)
	results = results.sort_values('isolate')
	return results


def generate_output(results, output_type):
	"""
	Summarise and prepare outputs
	"""
	if args.output_type == "tsv":
		print(results)
	elif args.output_type == 'summary':
		variant_sets = results.groupby('isolate')['Mutations']
		variant_sets = variant_sets.apply(list)
		for isolate, variant_set in variant_sets.iteritems():
			variant_set = sorted(variant_set,
           						 key=lambda x: (x.split(':')[1],
                                				x.split(':')[0]))
			print(f"{isolate}: {variant_set}")


if __name__ == '__main__':

	parser = argparse.ArgumentParser("Script to summarise ivar variants "
								     "in a convenient human readable manner")
	parser.add_argument('-i', '--input', required=True, nargs="+",
                        type=check_file,
						help="List of ivar_variants.tsv files to summarise")
	parser.add_argument('-t', '--output_type', choices=['tsv', 'summary'],
						 required=True, help="Output full concatenated "
                         "ivar variants file (tsv) or just a sorted summary "
                         "for each isolate (summary)")
	parser.add_argument('-r', '--reference_gff', default='data/MN908947_3.gff3',
                        type=check_file,
						 help="Path to MN908947.3 gff3 file")

	args = parser.parse_args()

	ref_db = gffutils.create_db(str(args.reference_gff), dbfn='ref.db',
                                force=True, keep_order=True,
                                merge_strategy='merge',
                                sort_attribute_values=True)

	results = parse_inputs(args.input, ref_db)

	generate_output(results, args.output_type)

