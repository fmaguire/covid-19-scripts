#!/bin/bash

#augur filter is really slow
python extract_seqs.py --nextmeta 2021_01_12/metadata_2021-01-11_16-49.tsv --nextfasta 2021_01_12/sequences_2021-01-11_09-53.fasta -q "pangolin_lineage=='B.1.1.28'" -o 2021_01_12_all_b1128 -r

mkdir -p augur

# remove short or seqs without month information
augur filter --sequences 2021_01_12_all_b1128_seqs.fasta --metadata 2021_01_12_all_b1128_metadata.tsv --exclude-ambiguous-dates-by month --min-length 2700 --output augur/filtered_seqs.fasta

# mask problematic sites and start/end from defaults in ncov
augur mask --mask-sites $(awk -F $'\t' '$1!~/^#/ {printf $2; printf " "}' < problematic_sites_sarsCov2.vcf) --sequences augur/filtered_seqs.fasta --mask-from-beginning 100 --mask-from-end 50 --output augur/masked_seqs.fasta

# align with mafft
augur align --sequences augur/masked_seqs.fasta --nthreads 8 --method mafft --reference-sequence ../ncov/defaults/reference_seq.fasta --fill-gaps --output augur/aligned.fasta

# build tree with iqtree matching ncov settings
augur tree --alignment augur/aligned.fasta --tree-builder-args '-ninit 10 -n 4' --output augur/raw_tree.nwk --nthreads 8

# refine tree matching ncov settings 
augur refine --timetree --tree augur/raw_tree.nwk --alignment augur/aligned.fasta --metadata 2021_01_12_all_b1128_metadata.tsv --output-tree augur/refined_tree.nwk --output-node-data augur/branch_lengths.json --root Wuhan/Hu-1/2019 Wuhan/WH01/2019 --clock-rate 0.0008 --clock-std-dev 0.0004 --coalescent "skyline" --date-inference marginal --divergence-unit mutations --date-confidence --no-covariance --clock-filter-iqd 4

# infer ancestral states and mutations
augur ancestral --tree augur/refined_tree.nwk --alignment augur/aligned.fasta --output-node-data augur/nt_muts.json --inference joint --infer-ambiguous

# translate nt to aa
augur translate --tree augur/refined_tree.nwk --ancestral-sequences augur/nt_muts.json --reference-sequence ../ncov/defaults/reference_seq.gb --output-node-data augur/aa_muts.json

# infer traits
augur traits --tree augur/refined_tree.nwk --metadata 2021_01_12_all_b1128_metadata.tsv --output augur/traits.json --columns "country_exposure" --confidence --sampling-bias-correction 2.5

# add clade info
augur clades --tree augur/refined_tree.nwk --mutations augur/nt_muts.json augur/aa_muts.json --clades ../ncov/defaults/clades.tsv --output-node-data augur/clades.json

augur frequencies --method kde --metadata 2022_01_12_all_b1128_metadata.tsv --tree augur/refined_tree.nwk --min-date 2020.0 --pivot-interval 1 --narrow-bandwidth 0.05 --proportion-wide 0.0 --minimal-frequency 0.01 --stiffness 20 --inertia 0.2 --output augur/tip-frequencies.json

augur frequencies --method diffusion --alignments augur/aligned.fasta --gene-names nuc --metadata 2021_01_12_all_b1128_metadata.tsv --tree augur/refined_tree.nwk --min-date 2020.0 --pivot-interval 1 --minimal-frequency 0.01 --stiffness 20 --inertia 0.2 --output augur/tip-frequencies.json

python ../ncov/scripts/assign-colors.py --ordering ../ncov/defaults/color_ordering.tsv --color-schemes ../ncov/defaults/color_schemes.tsv --output augur/colors.tsv --metadata 2021_01_12_all_b1128_metadata.tsv
python ../ncov/scripts/construct-recency-from-submission-date.py --metadata 2021_01_12_all_b1128_metadata.tsv --output augur/recency.json
python ../ncov/scripts/add_pangolin_lineages.py --tree augur/refined_tree.nwk --output augur/pangolin.json

augur export v2 --tree augur/refined_tree.nwk --metadata 2021_01_12_all_b1128_metadata.tsv --node-data augur/branch_lengths.json augur/nt_muts.json augur/aa_muts.json augur/clades.json augur/traits.json augur/recency.json --colors augur/colors.tsv --auspice-config ../ncov/my_profiles/b1128_lineage_build/my_auspice_config.json --lat-longs ../ncov/defaults/lat_longs.tsv --title "B.1.1.28 Build" --output augur/auspice.json --include-root-sequence


