#!/usr/bin/env python

from collections import defaultdict
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio import Data

## modified from https://github.com/jts/ncov-watch
class Variant:
    def __init__(self, contig, position, reference, alt):
        self.contig = contig
        self.position = position
        self.reference = reference
        self.alt = alt
        self.name = None
    def __repr__(self):
        return repr(",".join([str(x) for x in [self.contig, self.position,
                                              self.reference, self.alt]]))


## modified from type_variants script https://github.com/cov-ert/type_variants
def get_nuc_position_from_aa_description(CDS_dict, cds, aa_pos):
    """
    given a CDS (eg. S) and the number of an amino acid in it, get the
    1-based start position of that codon in Wuhan-Hu-1 ref coordinates
    nuc_pos is an integer which is 1-based start pos of codon
    """

    if cds.lower() not in CDS_dict.keys():
        sys.stderr.write("I don't know about cds: %s \n" % cds)
        sys.stderr.write("please use one of: %s" % ",".join(CDS_dict.keys()))
        sys.exit(1)

    #if cds.lower() == "orf1ab":
    #    if aa_pos <= 4401:
    #        parsed_cds = "orf1a"
    #    else:
    #        parsed_cds = "orf1b"
    #        aa_pos = aa_pos - 4401
    #else:
    #    parsed_cds = cds
    parsed_cds = cds

    cds_tuple = CDS_dict[parsed_cds.lower()]

    nuc_pos = cds_tuple[0] + ((aa_pos - 1) * 3)

    if nuc_pos > cds_tuple[1]:
        sys.stderr.write("invalid amino acid position for cds %s : %d" % (cds, aa_pos))
        sys.exit(1)

    return nuc_pos


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")



def convert_type_variants_to_vcf(reference: dict, type_variants_path: Path) -> list:
    """

    """

    converted_variants = []

    # the "backwards_table" function in biopython only returns one codon for
    # some reason so we have to make it ourselves
    codon_table = Data.CodonTable.standard_dna_table
    back_codon_table = defaultdict(list)
    for codon, aa in codon_table.forward_table.items():
        back_codon_table[aa].append(codon)

    print(back_codon_table)

    with open(type_variants_path) as fh:
        for variant in fh:
            variant = variant.strip().split(':')
            print(variant)
            if variant[0] == 'aa':
                orf = variant[1]
                ref_aa = variant[2][0]
                alt_aa = variant[2][-1]
                aa_pos = int(variant[2][1:-1])
                nuc_pos = get_nuc_position_from_aa_description(reference,
                                                               orf, aa_pos)

                ref_codon = reference['genome'][nuc_pos: nuc_pos + 3]

                # get all possible codons for alt aa
                alt_codons = back_codon_table[alt_aa]
                # for each one check what snps are required
                codon_positions = [nuc_pos, nuc_pos + 1, nuc_pos + 2]

                #
                for alt_codon in alt_codons:

                    for codon_pos, ref_codon_nt, alt_codon_nt in zip(codon_positions,
                                                                     ref_codon,
                                                                     alt_codon):
                        if ref_codon_nt != alt_codon_nt:
                            print(f"{ref_codon_nt}{codon_pos}{alt_codon_nt}")


                var = None
            elif variant[0] == 'del':
                # get the ref sequence for deletion length
                del_start = int(variant[1][1:]) - 1 * 3
                del_end = del_start + int(variant[2]) * 3
                ref_nt = reference['genome'][del_start - 1: del_end]
                alt_nt = ref_nt[0]
                var = Variant("MN908947.3", del_start, ref_nt, alt_nt)
            elif variant[0] == 'snp':
                ref_nt = variant[2][0]
                alt_nt = variant[2][-1]
                nt_pos = int(variant[2][1:-1])
                var = Variant("MN908947.3", nt_pos, ref_nt, alt_nt)
            else:
                raise ValueError("Invalid variant {variant}")

            converted_variants.append(var)

    return converted_variants

def convert_vcf_to_type_variants(vcf_path):
    """

    """
    pass


def parse_genbank(gbk_path: Path) -> dict:
    """
    Parse a gbk file for the reference genome into a dictionary for quick
    lookup
    """
    ref_data = {}
    for record in SeqIO.parse(gbk_path, 'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                gene = feature.qualifiers['gene'][0].lower()
                start = feature.location.start
                end = feature.location.end
                ref_data[gene] = (start, end)
        ref_data['genome'] = record.seq
    return ref_data


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Convert between type_variants format and "
                                     "vcf/ncov-watch format")

    parser.add_argument("-i", "--input", required=True, type=check_file,
                         help='Path to input vcf or type_variants config')
    parser.add_argument("-t", "--input_type", required=True,
                        choices=['vcf', 'type_variants_config'],
                        help="Is input a standard VCF (for ncov-watch) or a "
                             "a config file for type_variants")
    parser.add_argument("-g", "--ref_gbk", required=True, type=check_file,
                         help="Path to full genbank annotation for reference genome")

    args = parser.parse_args()

    reference = parse_genbank(args.ref_gbk)

    if args.input_type == 'type_variants_config':
        converted_vcf = convert_type_variants_to_vcf(reference, args.input)
        print(converted_vcf)

    elif args.input_type == "vcf":
        convert_type_variants = convert_vcf_to_type_variants(args.input,
                                                             reference)
        print(convert_type_variants)

