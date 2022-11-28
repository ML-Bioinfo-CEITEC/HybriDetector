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
- Run the app `snakemake --snakefile HybriDetector.smk --directory /path/to/the/directory/ --configfile configuration_file.json  --use-conda --conda-frontend mamba -p --res mem=100 -j 20`
`
