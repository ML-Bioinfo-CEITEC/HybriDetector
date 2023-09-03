# HybriDetector
HybriDetector is an automatic bioinformatics tool for the chimeric interactions detection between variety of noncoding RNAs and their genomic targets. HybriDetector allows to obtain hybrid interactions from binding NGS experiments namely CLASH (crosslinking, ligation and sequencing of hybrids), CLEAR-CLIP (covalent ligation of endogenous Argonaute-bound RNAs) or new miR-eCLIP and is also suitable to analyze simple RNA-RNA binding sequencing data based on regular CLIP-seq protocols and other protocols, where intermolecular ligation might naturally occur during the library preparation steps. HybriDetector deals with all steps of data analysis starting with preprocessed trimmed fastqc files and producing the table containing obtained molecular interactions between noncoding RNAs and their genomic targets together with frequency, annotations and predicted hybrid secondary structure folding.     

## Method overview

Sequencing reads resulting from experiments mentioned above consists of two distinct interacting RNA molecules, partially digested and connected by intermolecular ligation at one end. At both ends of molecular construct the adaptors necessary for high throughput sequencing might be present. To be able to conclusively determine PCR duplicates introduced during the library preparation for sequencing, the 5’ end adapter are usually enriched with Unique Molecular Identifiers (UMI). To obtain and annotate highly accurate chimeric reads within the provided data `HybriDetector` needs to distinguish between a huge amount of very similar sequences of distinct non-coding RNAs. This makes it very sensitive to every single base difference in the input reads. Therefore to produce more accurate  very precisly trimmed and cleaned input fastq files are required   

HybriDetector bioinformatics solution begins with the alignment of input fastq files to the whole human genome using very strict alignment settings to identify single genomic reads. Mapping to the whole genome gives us an outstanding opportunity to obtain all non common and purely described interactions within introns or even in intergenic regions. All reads uniquely aligned to the human genome are separated and considered as *non-chimeric single genomic* reads apart from those originating from our masked regions. Masked region consists of all non-coding RNAs belonging to one of the following biotypes: miRNA, tRNA, rRNA, snoRNA, vaultRNA and YRNA regions gathrered from NCBI, Ensembl, UCSC, RepeatMasker, GtRNAdb, SILVA and DASHR databases.

All of the unmapped reads together with reads from masked regions subject to a another round of the alignment to every single non-coding RNA list mentioned before, serving to determine both the single non-coding reads as well as chimeric candidates. Reads fully covering any non-coding RNA without any significant sequence overhang are separated and considered as single non-coding RNA non-chimeric reads. Reads partially covering any non-coding RNA and at the same time showing significant sequence overhang are moved further and considered as potential chimeric candidates.

Unmapped overhang parts of reads obtained after mapping to individual non-coding RNA databases are further mapped to the whole genome. Perfectly mapped overhang sequence reads determine the genomic part of potential chimeric reads.

After the chimeras are called, all duplicated records are removed using UMI’s and postprocess filters i.e. no mismatch in alignments, support of chimeric interaction by single reads alignment, alignment lengths of individual chimeric read parts, etc. are applied to further refine called chimeras. To reduce the redundancy of called chimeras we collapsed all chimeras belonging to the specific non-coding RNA and targeting the same genomic region with some minor discrepancies into one chimera.

Next to the filters, interactions and base pairing within every single called chimera is predicted by folding the chimeric sequence with `VienaRNA`. At the same time genomic features, RepeatMasker annotation and number of unique and duplicated reads supporting individual chimera are assigned to every obtained hybrid.        


## Installation
To install HybriDetector, you need the following prerequisites:
- [Anaconda](https://www.anaconda.com/products/individual) - Python package manager for safe installation 
- UNIX based OS
- Minimum of 6 CPU cores and 24GB RAM to analyze ~ 35M of 120 bp long sequenced reads

### Manual installation
If you do not have [Anaconda](https://www.anaconda.com/distribution/) installed on your computer, please do so first. 
- Download the latest release from [the repository](https://github.com/ML-Bioinfo-CEITEC/CLASH) with `git clone https://github.com/ML-Bioinfo-CEITEC/HybriDetector.git`
- Go to the project directory `cd HybriDetector`
- Recreate environment from yml file `conda env create -f HybriDetector_env.yml`
- Activate the environment `conda activate HybriDetector_env`
- Create the folder for the input file/s called "data" `mkdir data`
- Run the app `snakemake --snakefile HybriDetector.smk --configfile configuration_file.json  --use-conda --conda-frontend mamba -p --res mem=100 -j 20`
- All of the references, databases and annotations are either directly colned with the repository or downloaded during the first run of the tool.

## Hybrid detection
Required input is already preprocessed zipped fastq file with precisely removed sequencíng barcodes and adapters and with suffix ".fastq.gz". In case the UMIs (Unique molecular identifiers) are used in a sequencíng library, should be already extracted and moved from the sequence to the read header. HybriDetector is capable to analyze single sample as well as multiple replicates of particular samples within one run. HybriDetector can currently be used for the human genome only.

### Arguments 

`--input_sample` - Fastq file names here written without the ".fastq.gz" suffix separated by comma.

`--read_length` - [75] Maximal length of sequenced read within input fastq file (integer).

`--is_umi` - [FALSE] Specify whether your library contains extracted UMIs in the read header (boolean).

`--map_perc_single_genomic` - [0.85] Specify the percentage when the alignment of single genomic non-chimeric reads will be reported. Only if the ratio of "aligned length / read length" is higher than or equal to this
value, the alignemnt will be output.

`--map_perc_softclip` - [0.75] Specify the percentage when the alignment of genomic part of the chimeric reads will be reported. Only if the ratio of "aligned length / read length" is higher than or equal to this
value, the alignemnt will be output.

`--cores - [6] Number of provided CPUs.`

`--ram - [24] Number of provided RAM in GB.`

### Usage examples

`python HybriDetector.py --input_sample SAMPLE1.fastq.gz,SAMPLE2.fastq.gz --read_length 75 --is_umi FALSE --map_perc_single_genomic 0.85 --map_perc_softclip 0.75 --cores 6 --ram 24` 

Command above starts analysis based on two already trimmed replicates where there are UMIs not used requiring 6 CPU cores and 24 GB of RAM.

### Output files

All of the outputs are stored within a folder `HybriDetector/hyb_pairs/`.

Two main branches of the output file are created based on the alignment uniqueness of the non-coding RNA part of the detected chimeric interaction. In a case, where the identified non-coding RNA part of the hybrid uniquely align to just one non-coding RNA, these interactions are treated as `unique`. All the rest are reported as `ambiguous`. All of the hybrid outputs are named considering that. The list of important outputs and their description:

#### `{sample}.unified_length_*_types_*_high_confidence_finalout.tsv`
The main and final output containing all identified hybrid interactions together with their annotations. Description of individual columns:

`ID` - unique identifier of the obtained hybrid interaction.

`Genomic fragment sequence` - Sequence of genomic part of the hybrid interaction anchored to the middle and prolong equally at both ends to the length of 50 nt (in case obtained sequence was longer than 50 nt, it is cut equally from both sides to the lenght of 50 nt). 

`Driver fragment sequence` - Sequence of non-coding RNA part of the hybrid interaction.

`Target first in chimera` - True when the order of the hzbrid parts was target-driver (mRNA-noncodingRNA).

`N chimeric reads` - Number of reads supporting obtained hybrid interaction.

`N deduplicated chimeric reads` - Number of reads supporting obtained hybrid interaction after deduplication with UMIs (available only when used with UMIs).

`Chimera in replicates` - List of replicates, where particular hybrid interactoin occured (available only when replicates provided).

`N non-chimeric reads` - Number of *non-chimeric single genomic* reads supporting obtained hybrid interaction.

`Coverage non-chimeric reads` - Normalized coverage of *non-chimeric single genomic* reads supporting obtained hybrid interaction.

`Chromosome` - Chromosome where hybrid interaction occured.

`Start` - Start of the loci, where genomic part of the hybrid was mapped to.

`End` - End of the loci, where genomic part of the hybrid was mapped to.

`Strand` - Whether the obtained hybrid interaction occured on the forward or reverse genomic strand.

`Overlapping gene name` - Gene name in which hybrid interaction was detected.

`Overlapping gene biotype` - Biotype of the gene (protein coding, antisense, etc.) in which hybrid interaction was detected.

`Overlapping gene feature` - At which genomic feature (exon, intron, etc.) hybrid interaction was detected.

`Overlapping repeatmasker` - Whether the loci, where hybrid interaction was detected overlap with any record in RepeatMasker annotation.

`Overlapping repeatmasker family` - Into which family the RepeatMasker record belongs to (in case the RepeatMasker overlap is detected).

`Driver name` - Name of mapped non-coding RNA.

`Driver type` - Type of the non-coding RNA (from the list above) obtained within hybrid interaction.

`miRNA Family` - Family to which microRNA belongs to (in case that non-coding RNA type is microRNA).

`Driver alignment no mismatch` - Whether the gathered sequence of non-coding hybrid part mapped perfectly to its reference.

`Cofold structure` - Folded secondary structure of obtained hybrid interaction predicted utilizing ViennaRNA toolkit.

`Cofold MFE` - Free energyof predicted folded secondary structure.

#### `{sample}.unified_length_*_types_*_finalout.tsv`
This output table contain the same columns as the main output file, the difference rests in not using of filtering based on "Driver alignment no mismatch" column. Within this output file, there are both these hybrid interactions where seqeunce of non-coding RNA part mapped perfectly to the reference and these where aligment shows up to 2 mismatches.  

#### `{sample}.all_hyb_pairs.tsv`
First and raw intermediate output. Contains all af the obtained hybrid interaction pairs, without any filtering and deduplication. Description of individual columns:

`read_name` - Original name of sequenced read, where hybrid interaction was identified.

`flag.g` - Bitwise FLAG assigned by STAR aligner for the mapping of genomic part of the hybrid interaction. 

`chr.g` - Chromosome where mapping of genomic part of the hybrid interaction occured.

`pos.g` - Start of the loci, where genomic part of the hybrid was mapped to.

`qual.g` - Mapping quality assigned by STAR aligner for the mapping of genomic part of the hybrid interaction. 

`cigar.g` - CIGAR string assigned by STAR aligner (describes the precision of mapping) for the mapping of genomic part of the hybrid interaction.

`seq.g` - Sequence of genomic part of the hybrid interaction.

`gene_name` - Gene name in which hybrid interaction was detected.

`strand` - Whether the obtained hybrid interaction occured on the forward or reverse genomic strand.

`smallrna` - Name of mapped non-coding RNA.

`pos.m` - Start of the loci, where non-coding RNA part of the hybrid was mapped to.

`qual.m` -  Mapping quality assigned by STAR aligner for the mapping of non-coding RNA part of the hybrid interaction.

`cigar.m` - CIGAR string assigned by STAR aligner (describes the precision of mapping) for the mapping of non-coding RNA part of the hybrid interaction.

`seq.m` - Sequence of non-coding RNA part of the hybrid interaction. 

`forward_dir` - True when the order of the hybrid parts was target-driver (mRNA-noncodingRNA).

`UMI` - Sequence of UMI for particular sequenced read used further for the deduplication.

#### Folders in a path `HybriDetector/hyb_pairs/cofold*`
Within these output folders can be found folded secondary structure predicted by utilizing of ViennaRNA toolkit in a postscript format. Each predicted secondary structure can be matched with obtained hybrid interaction by concatenation of the `target` and `driver` sequence with `&` as separator, as all predicted secondary structure are named in such format.

### Tool list
Tools incorporated into the HybriDetector pipeline necesarry to gather the hybrid interactions are `FastQC, Samtools, STAR, deepTools, UMI-tools, Picard, r-data.table, r-seqinr,r-stringi, r-proxy, bioconductor-biostrings, bioconductor-rtracklayer, bioconductor-genomicranges, bioconductor-bsgenome.hsapiens.ucsc.hg38, bioconductor-mirbaseconverter`. 