# HybriDetector
HybriDetector is an automatic bioinformatics tool for the chimeric interactions detection between variety of noncoding RNAs and their genomic targets. HybriDetector allows to obtain hybrid interactions from binding NGS experiments namely CLASH (crosslinking, ligation and sequencing of hybrids), CLEAR-CLIP (covalent ligation of endogenous Argonaute-bound RNAs) or new chimeric eCLIP and is also suitable to analyze simple RNA-RNA binding sequencing data based on regular CLIP-seq protocols and other protocols, where intermolecular ligation might naturally occur during the library preparation steps. HybriDetector deals with all steps of data analysis starting with preprocessed trimmed fastqc files and producing the table containing obtained molecular interactions between noncoding RNAs and their genomic targets together with frequency, annotations and predicted hybrid secondary structure folding.     

## Installation
To install HybriDetector, you need the following prerequisites:
- [Anaconda](https://www.anaconda.com/products/individual) - Python package manager for safe installation 
- UNIX based OS
- Minimum of 6 CPU cores and 24GB RAM to analyze ~ 35M of 120 bp long sequenced reads

### Manual installation
If you do not have [Anaconda](https://www.anaconda.com/distribution/) installed on your computer, please do so first. 
- Download the latest release from [the repository](https://github.com/ML-Bioinfo-CEITEC/CLASH) with `git clone https://github.com/ML-Bioinfo-CEITEC/HybriDetector.git`
- Go to the project directory `cd HybriDetector`
- Recreate environment from yml file `conda env create -f snakemake.yml`
- Activate the environment `conda activate snakemake`
- Create the folder for the input file/s called "data" `mkdir data`
- Run the app `snakemake --snakefile HybriDetector.smk --directory /path/to/the/directory/ --configfile configuration_file.json  --use-conda --conda-frontend mamba -p --res mem=100 -j 20`

## Hybrid detection
Required input is already preprocessed zipped fastq file with precisely removed sequencíng barcodes and adapters and with suffix ".fastq.gz". In case the UMIs (Unique molecular identifiers) are used in a sequencíng library, should be already extracted and moved from the sequence to the read header. HybriDetector is capable to analyze single sample as well as multiple replicates of particular samples within one run. 

### Arguments 

`--sample` - fastq file names here written without the ".fastq.gz" suffix sepparated by space
`--read_length` - [75] maximal lenght of seqeunced read within input fastq file
`--is_umi` - [FALSE] specify whether you library contains extracted UMIs in the read header
`--map_perc_single_genomic`
`--map_score_single_genomic`
`--map_perc_softclip`
`--map_score_softclip`