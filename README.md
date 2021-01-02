# SARS-CoV-2 analysis scripts

Repository to keep track of miscellaneous SARS-CoV-2 analysis scripts

## Summarise Ivar

Quick script to translate and provide quick human readable breakdown of 
mutations from ivar variants.  Hopefully brings some of the convenience of
the nextclade summaries using the more robust ivar inferences instead of
relying on minimap2 alignments.

Very loosely based on parts of [type_variants](github.com/cov-ert/type_variants).

### Installation

Requires gffutils and pandas to work:

    conda create -n summarise_ivar gffutils pandas
    conda activate summarise_ivar

### Usage

    python summarise_ivar.py --input test/test1_ivar_variants.tsv test/test2_ivar_variants.tsv --output_type summary
    >>> test/test1_ivar_variants.tsv: ['aa:N:R203K', 'aa:N:G204R', 'snp:N:G28882A', 'aa:S:D614G', 'snp:non-coding:C241T', 'aa:orf1ab:P4715L', 'snp:orf1ab:C3037T', 'snp:orf1ab:T7288C', 'snp:orf1ab:A14199G']
    >>> test/test2_ivar_variants.tsv: ['aa:N:R203K', 'aa:N:G204R', 'snp:non-coding:C241T', 'snp:orf1ab:C3037T', 'snp:orf1ab:T7288C']


## Lineage Assignments

Generate pangolin and nextclade lineage assignments before collating results
into a single table containing version metadata for tools and data dependencies

### Installation

Requires nextclade and pangolin to work:

    conda env create -f https://raw.githubusercontent.com/cov-lineages/pangolin/master/environment.yml
    conda activate pangolin
    pip install git+https://github.com/cov-lineages/pangolin
    pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade
    pip install git+https://github.com/cov-lineages/lineages.git --upgrade
    conda install -y nodejs 
    npm install --global @neherlab/nextclade

### Usage 

    python assign_lineages.py --input_genomes test/genomes.fa --output lineages.tsv

