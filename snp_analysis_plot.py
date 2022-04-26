#!/usr/bin/env python
import pandas as pd
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import vcf
sns.set_style('whitegrid')
sns.set_palette('colorblind')


def parse_vcf_snpeff_record(record, sample):
    parsed_records = []
    for ix, alt in enumerate(record.ALT):
        parsed_alt = {}

        # get var only for that alt
        var_list = [ann for ann in record.INFO['ANN'] if ann.split('|')[0] == alt]
        #remove up/downstream/intergenic mutations
        var_list = [ann for ann in var_list if ann.split('|')[1] not in ['upstream_gene_variant',
                                                                                 'downstream_gene_variant']]
        if len(var_list) > 1:
            print("var_list has multiple changes", var_list)
            assert False
        elif len(var_list) == 0:
            print("var_list has no annotation", var_list)
            assert False
        else:
            ann = var_list[0].split('|')
            parsed_alt['Mutation Effect'] = ann[1]
            parsed_alt['Mutation Gene'] = ann[3]
            parsed_alt['Nucleotide Mutation'] = ann[9]

            if ann[1] == 'intergenic_region':
                parsed_alt['Protein Mutation'] = f"No Protein Effect ({ann[9]})"
            elif ann[1] == 'synonymous_variant':
                parsed_alt['Protein Mutation'] = ann[3] + ": synonymous " + ann[10]
            else:
                parsed_alt['Protein Mutation'] = ann[3] + ":" + ann[10]

        parsed_alt['Sample'] = sample
        parsed_alt['Genome Position'] = record.POS
        parsed_alt['Allele Read Count'] = record.INFO['AO'][ix]
        parsed_alt['Total Read Count'] = record.INFO['DP']
        parsed_alt['% Reads Supporting Allele'] = record.INFO['VAF'][ix] * 100

    parsed_records.append(parsed_alt)
    return parsed_records

def plot_allele_pres_absence(variants_subset, title, savepath, all_mutations=False):

    coverage_thresold = 50

    variant_order = variants_subset.sort_values("Genome Position")['Protein Mutation'].unique()
    sample_order = variants_subset['Sample'].unique()
    variant_presence_absence = pd.crosstab(variants_subset['Sample'], variants_subset['Protein Mutation'])
    variant_presence_absence = variant_presence_absence.loc[sample_order, variant_order].T

    variant_percentage_alleles = variants_subset[variants_subset['Protein Mutation'].isin(variant_presence_absence.index)]
    variant_percentage_alleles = variant_percentage_alleles.pivot('Sample', 'Protein Mutation', '% Reads Supporting Allele').fillna(0)
    variant_order = [x for x in variants_subset.sort_values('Genome Position')['Protein Mutation'].unique() if x in variant_percentage_alleles.columns]
    variant_percentage_alleles = variant_percentage_alleles.loc[sample_order, variant_order].T


    # drop any all >99%
    if not all_mutations:
        all_high_coverage_variants = [var for var in variant_percentage_alleles.T if (variant_percentage_alleles.T[var] > 95).all() ]
        variant_percentage_alleles = variant_percentage_alleles.drop(all_high_coverage_variants)

    # low coverage mutations (possible dropout)
    low_coverage = variants_subset[variants_subset['Allele Read Count'] < coverage_thresold]['Protein Mutation'].unique()

    #variant_percentage_alleles = variant_percentage_alleles.rename(index=possible_dropout)

    #return variant_percentage_alleles
    if len(variant_percentage_alleles) < 10:
        fig, ax = plt.subplots(figsize=(6,8))
    elif len(variant_percentage_alleles) > 30:
        fig, ax = plt.subplots(figsize=(6,16))
    else:
        fig, ax = plt.subplots(figsize=(6,12))



    #variant_percentage_alleles = variant_percentage_alleles[sorted(variant_percentage_alleles.columns, key=lambda x: x.split('-')[1])]

    ax.set_title(title)
    sns.heatmap(variant_percentage_alleles, vmin=0, vmax=100, linewidths=.1, ax=ax, xticklabels=True, cmap="mako_r", yticklabels=True,
                cbar_kws={'label': '% Reads Supporting Allele'})
    ax.set_ylabel(f"Protein Mutation\n(Any Sample <{coverage_thresold}X Coverage in Red)")
    ax.set_xlabel(f"Samples\n(Original Genomes in Bold)")

    # colour OG and low coverage labels
    for label in ax.get_yticklabels():
        if label.get_text() in low_coverage:
            label.set_color('red')

     # colour OG and low coverage labels
    for label in ax.get_xticklabels():
        if "OG" in label.get_text():
            label.set_fontweight('bold')
    #return label
    #return fig, ax
    print(f"Writing figure to {savepath}")
    plt.savefig(savepath, dpi=300, bbox_inches='tight', facecolor='white', transparent=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Tool to summarise SNPs across samples")
    parser.add_argument("-n", "--name", required=True, help="Analysis name for output/titles")
    parser.add_argument("vcfs", nargs="+", help="List of SnpEff annotated VCFs")
    args = parser.parse_args()

    vcfs = []
    for vcf_fp in args.vcfs:
        vcf_fp = Path(vcf_fp)
        if not vcf_fp.exists():
            raise ValueError(f"{vcf_fp} does not exist")
        else:
            vcfs.append((vcf_fp.name.replace('.vcf', '').replace('.ann', ''), vcf_fp.resolve()))

    parsed_records = []
    for sample, vcf_fp in vcfs:
        with open(vcf_fp) as fh:
            for record in vcf.Reader(fh):
                parsed_records.extend(parse_vcf_snpeff_record(record, sample))

    variants = pd.DataFrame(parsed_records)
    plot_allele_pres_absence(variants, args.name, f'{args.name}_passaging_all.png', all_mutations=True)


