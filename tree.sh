#!/bin/bash

#augur filter is really slow
python extract_seqs.py --nextmeta 2021_01_12/metadata_2021-01-11_16-49.tsv --nextfasta 2021_01_12/sequences_2021-01-11_09-53.fasta -q "pangolin_lineage=='B.1.1.28'" -o 2021_01_12_all_b1128 -r

mkdir -p augur

# remove short or seqs without month information
augur filter --sequences 2021_01_12_all_b1128_seqs.fasta --metadata 2021_01_12_all_b1128_metadata.tsv --exclude-ambiguous-dates-by month --min-length 2700 --output augur/filtered_seqs.fasta

# mask problematic sites and start/end from defaults in ncov
augur mask --mask-sites $(awk -F $'\t' '$1!~/^#/ {printf $2; printf " "}' < problematic_sites_sarsCov2.vcf) --sequences augur/filtered_seqs.fasta --mask-from-beginning 100 --mask-from-end 50 --output augur/masked_seqs.fasta

# align with mafft
augur align --sequences augur/masked_seqs.fasta --nthreads 8 --method mafft --reference-sequence data/MN908947.fa --fill-gaps --output augur/aligned.fasta

# build tree with iqtree matching ncov settings
augur tree --alignment augur/aligned.fasta --tree-builder-args '-ninit 10 -n 4' --output augur/raw_tree.nwk --nthreads 8

# refine tree matching ncov settings 
augur refine --tree augur/raw_tree.nwk --alignment augur/aligned.fasta --metadata 2021_01_12_all_b1128_metadata.tsv --output-tree augur/refined_tree.nwk --output-node-data augur/branch_lengths.json --root Wuhan/Hu-1/2019 Wuhan/WH01/2019 --clock-rate 0.0008 --clock-std-dev 0.0004 --coalescent "skyline" --date-inference marginal --divergence-unit mutations --date-confidence --no-covariance --clock-filter-iqd 4

# infer ancestral states and mutations
augur ancestral --tree augur/refined_tree.nwk --alignment augur/aligned.fasta --output-node-data augur/nt_muts.json --inference joint --infer-ambiguous

# translate nt to aa
augur translate --tree augur/refined_tree.nwk --ancestral-sequences augur/nt_muts.json --reference-sequence data/reference_seq.gb --output-node-data augur/aa_muts.json

# infer traits
augur traits --tree augur/refined_tree.nwk --metadata 2021_01_12_all_b1128_metadata.tsv --output augur/traits.json --columns "country_exposure" --confidence --sampling-bias-correction 2.5

# add clade info
augur clades --tree augur/refined_tree.nwk --mutations augur/nt_muts.json augur/aa_muts.json --clades data/clades.tsv --output-node-data augur/clades.json

#


