import math
import subprocess
import json
import re
import os.path
import glob
import pandas as pd
from snakemake.utils import R
from snakemake.utils import report
from os.path import split

###################
## CLASH analysis ##
###################

## REFERENCES AND BEDS
cfg = pd.DataFrame(config)

TMP_DIR = "/mnt/ssd/ssd_1/tmp/"
#REF_DIR = define_variable(cfg, "REF_DIR")
REF = "/DBs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REF_MIR = "/DBs/mirna/hsa_mature.fa"
REF_TRNA = "/DBs/tRNA_final.fa"
REF_RRNA = "/DBs/rrna/rRNA_ensembl_dedup_final.fa"
REF_SNORNA = "/DBs/snorna/snorna_dashr_dedup_final.fa"
REF_TRNA_ENS = "/DBs/trna/ensembl/tRNA_ensembl_dedup_final.fa"
REF_TRNA_GTRNADB = "/DBs/trna/gtrnadb/tRNA_gtrnadb_dedup_final.fa"
REF_TRNA_NCBI = "/DBs/trna/ncbi/tRNA_ncbi_dedup_final.fa"
REF_TRNA_UCSC= "/DBs/trna/ucsc/tRNA_ucsc_dedup_final.fa"
REF_VAULTRNA = "/DBs/vaultrna/vaultRNA_dedup_final.fa"
REF_YRNA = "/DBs/YRNA_final.fa"
REF_YRNA_ENS = "/DBs/yrna/ensembl/YRNA_ensembl_dedup_final.fa"
REF_YRNA_UCSC = "/DBs/yrna/ucsc/YRNA_ucsc_dedup_final.fa"
NON_CODING_RNA_BED = "/non_coding_rna_sorted_merged.bed"
REPEATMASKER = "/DBs/repeatmasker.bed"
TRANSCRIPTS = "/DBs/GRCh38-p10_transcripts_of_longest_transcripts.gtf"
FEATURES = "/DBs/GRCh38-p10_longest_transcripts_subtracted_final.gtf" 


DATA_DIR = "/mnt/ssd/ssd_3/temp/vasek/CLASH"
SAMPLES = cfg["Sample"]
print(len(SAMPLES))
wildcard_constraints:
    type="(unique|ambiguous)"



def all_inputs(wildcards):
    if len(SAMPLES) > 1:
        return expand(DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_unified_length_all_types_{type}_high_confidence.tsv", sample = "Merged", type = ["unique","ambiguous"])
    else: 
        return expand(DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_unified_length_all_types_{type}_high_confidence.tsv", sample = SAMPLES, type = ["unique","ambiguous"])
rule all:
    input:
        pair_collapsed_unified = all_inputs
        #pair_collapsed_unified = expand(DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_unique.tsv", sample = SAMPLES)

rule fastqc:
	input:	
		reads = DATA_DIR + "/preprocessed/{sample}.fastq.gz"
	output:
		html = DATA_DIR + "/raw_fastq_qc/{sample}.fastqc.html"
	log:
		run = DATA_DIR + "/sample_logs/{sample}/raw_fastq_qc.log"
	params:	
		prefix = DATA_DIR + "/raw_fastq_qc/",
		html = DATA_DIR + "/raw_fastq_qc/{sample}_fastqc.html"
	threads:	1
	conda:	"wraps/fastqc/env.yaml"
	script:	"wraps/fastqc/script.py"


rule STAR_gen_index:
    input:  gen = DATA_DIR + REF,       
    output: SAindex = DATA_DIR + "/index/STAR/SAindex",
    params: dir = DATA_DIR + "/index/STAR",
            log = DATA_DIR + "/index/STAR/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR.indexation_run.log",
    threads:    40
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_mirna_gen_index:
    input:  gen = DATA_DIR + REF_MIR,       
    output: SAindex = DATA_DIR + "/index/STAR_mirna/SAindex",
    params: dir = DATA_DIR + "/index/STAR_mirna",
            log = DATA_DIR + "/index/STAR_mirna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_mirna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_rrna_gen_index:
    input:  gen = DATA_DIR + REF_RRNA,       
    output: SAindex = DATA_DIR + "/index/STAR_rrna/SAindex",
    params: dir = DATA_DIR + "/index/STAR_rrna",
            log = DATA_DIR + "/index/STAR_rrna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_rrna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_snorna_gen_index:
    input:  gen = DATA_DIR + REF_SNORNA,       
    output: SAindex = DATA_DIR + "/index/STAR_snorna/SAindex",
    params: dir = DATA_DIR + "/index/STAR_snorna",
            log = DATA_DIR + "/index/STAR_snorna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_snorna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ens_gen_index:
    input:  gen = DATA_DIR + REF_TRNA_ENS,       
    output: SAindex = DATA_DIR + "/index/STAR_trna_ens/SAindex",
    params: dir = DATA_DIR + "/index/STAR_trna_ens",
            log = DATA_DIR + "/index/STAR_trna_ens/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ens.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_gtrnadb_gen_index:
    input:  gen = DATA_DIR + REF_TRNA_GTRNADB,       
    output: SAindex = DATA_DIR + "/index/STAR_trna_gtrnadb/SAindex",
    params: dir = DATA_DIR + "/index/STAR_trna_gtrnadb",
            log = DATA_DIR + "/index/STAR_trna_gtrnadb/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_gtrnadb.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ncbi_gen_index:
    input:  gen = DATA_DIR + REF_TRNA_NCBI,       
    output: SAindex = DATA_DIR + "/index/STAR_trna_ncbi/SAindex",
    params: dir = DATA_DIR + "/index/STAR_trna_ncbi",
            log = DATA_DIR + "/index/STAR_trna_ncbi/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ncbi.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ucsc_gen_index:
    input:  gen = DATA_DIR + REF_TRNA_UCSC,       
    output: SAindex = DATA_DIR + "/index/STAR_trna_ucsc/SAindex",
    params: dir = DATA_DIR + "/index/STAR_trna_ucsc",
            log = DATA_DIR + "/index/STAR_trna_ucsc/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ucsc.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_vaultrna_gen_index:
    input:  gen = DATA_DIR + REF_VAULTRNA,       
    output: SAindex = DATA_DIR + "/index/STAR_vaultrna/SAindex",
    params: dir = DATA_DIR + "/index/STAR_vaultrna",
            log = DATA_DIR + "/index/STAR_vaultrna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_vaultrna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_yrna_ens_index:
    input:  gen = DATA_DIR + REF_YRNA_ENS,       
    output: SAindex = DATA_DIR + "/index/STAR_yrna_ens/SAindex",
    params: dir = DATA_DIR + "/index/STAR_yrna_ens",
            log = DATA_DIR + "/index/STAR_yrna_ens/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_yrna_ens.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_yrna_ucsc_gen_index:
    input:  gen = DATA_DIR + REF_YRNA_UCSC,       
    output: SAindex = DATA_DIR + "/index/STAR_yrna_ucsc/SAindex",
    params: dir = DATA_DIR + "/index/STAR_yrna_ucsc",
            log = DATA_DIR + "/index/STAR_yrna_ucsc/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_yrna_ucsc.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule alignment_single_genomic:
    input:  fastq = DATA_DIR + "/preprocessed/{sample}.fastq.gz",
            genome = DATA_DIR + REF,
            index = DATA_DIR + "/index/STAR/SAindex",
    output: bam = DATA_DIR + "/mapped/to_genome/{sample}.to_genome.bam",
            bai = DATA_DIR + "/mapped/to_genome/{sample}.to_genome.bam.bai",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment.log"
    threads: 20
    resources:  mem = 34
    params: prefix = DATA_DIR + "/mapped/to_genome/{sample}/{sample}.",
            read_len = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min(), # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html
            map_perc = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_perc_single_genomic"].min(),
            map_score = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_score_single_genomic"].min(),
            snakemake_dir = workflow.basedir + "/../",
            tmpd = TMP_DIR,
    conda:  "wraps/alignment/env.yaml"
    script: "wraps/alignment/script.py"

rule filter_single_genomic:
    input:  bam = DATA_DIR + "/mapped/to_genome/{sample}.to_genome.bam",
            bai = DATA_DIR + "/mapped/to_genome/{sample}.to_genome.bam.bai",
            bed = DATA_DIR + NON_CODING_RNA_BED
    output: bam = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic.bam",
            bam_non_genomic = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.bam",
            fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz"
    log:    run = DATA_DIR + "/sample_logs/{sample}/filter_single_genomic.log"
    threads: 20
    resources:  mem = 34
    params: prefix = DATA_DIR + "/mapped/genomic/{sample}.",
            prefix_non_genomic = DATA_DIR + "/mapped/non_genomic/{sample}.",
            
    conda:  "wraps/filter_single_genomic/env.yaml"
    script: "wraps/filter_single_genomic/script.py"

rule deduplicate_single_genomic:
    input: bam = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic.bam"
    output: bam = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic_dedup.bam",
            bed = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm.bedgraph"
    log:    run = DATA_DIR + "/sample_logs/{sample}/deduplicate_single_genomic.log"
    threads: 20
    resources:  mem = 34
    conda:  "wraps/deduplicate_single_genomic/env.yaml"
    script: "wraps/deduplicate_single_genomic/script.py"

rule alignment_mirna:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_MIR,
            index = DATA_DIR + "/index/STAR_mirna/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/mirna_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_mirna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/mirna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_rrna:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_RRNA,
            index = DATA_DIR + "/index/STAR_rrna/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/rrna_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_rrna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/rrna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_snorna:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_SNORNA,
            index = DATA_DIR + "/index/STAR_snorna/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/snorna_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_snorna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/snorna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ens:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_TRNA_ENS,
            index = DATA_DIR + "/index/STAR_trna_ens/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/trna_ens_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_trna_ens.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/trna_ens_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_gtrnadb:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_TRNA_GTRNADB,
            index = DATA_DIR + "/index/STAR_trna_gtrnadb/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/trna_gtrnadb_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_trna_gtrnadb.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/trna_gtrnadb_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ncbi:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_TRNA_NCBI,
            index = DATA_DIR + "/index/STAR_trna_ncbi/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/trna_ncbi_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_trna_ncbi.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/trna_ncbi_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ucsc:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_TRNA_UCSC,
            index = DATA_DIR + "/index/STAR_trna_ucsc/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/trna_ucsc_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_trna_ucsc.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/trna_ucsc_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_vaultrna:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_VAULTRNA,
            index = DATA_DIR + "/index/STAR_vaultrna/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/vaultrna_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_vaultrna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/vaultrna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_yrna_ens:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_YRNA_ENS,
            index = DATA_DIR + "/index/STAR_yrna_ens/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/yrna_ens_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_yrna_ens.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/yrna_ens_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_yrna_ucsc:
    input:  fastq = DATA_DIR + "/mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = DATA_DIR + REF_YRNA_UCSC,
            index = DATA_DIR + "/index/STAR_yrna_ucsc/SAindex",
    output: bam = DATA_DIR + "/mapped/non_genomic/yrna_ucsc_alignment/{sample}.mapped.filtered.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_yrna_ucsc.log"
    threads: 20
    resources:  mem = 50
    params: prefix = DATA_DIR + "/mapped/non_genomic/yrna_ucsc_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule disambiguate_and_filter:
    input:  mirna_sam = DATA_DIR + "/mapped/non_genomic/mirna_alignment/{sample}.mapped.filtered.sam",
            rrna_sam = DATA_DIR + "/mapped/non_genomic/rrna_alignment/{sample}.mapped.filtered.sam",
            snorna_sam = DATA_DIR + "/mapped/non_genomic/snorna_alignment/{sample}.mapped.filtered.sam",
            trna_ens_sam = DATA_DIR + "/mapped/non_genomic/trna_ens_alignment/{sample}.mapped.filtered.sam",
            trna_gtrnadb_sam = DATA_DIR + "/mapped/non_genomic/trna_gtrnadb_alignment/{sample}.mapped.filtered.sam",
            trna_ncbi_sam = DATA_DIR + "/mapped/non_genomic/trna_ncbi_alignment/{sample}.mapped.filtered.sam",
            trna_ucsc_sam = DATA_DIR + "/mapped/non_genomic/trna_ucsc_alignment/{sample}.mapped.filtered.sam",
            vaultrna_sam = DATA_DIR + "/mapped/non_genomic/vaultrna_alignment/{sample}.mapped.filtered.sam",
            yrna_ens_sam = DATA_DIR + "/mapped/non_genomic/yrna_ens_alignment/{sample}.mapped.filtered.sam",
            yrna_ucsc_sam = DATA_DIR + "/mapped/non_genomic/yrna_ucsc_alignment/{sample}.mapped.filtered.sam",
    output: bam = DATA_DIR + "/mapped/noncoding_rna_chim/{sample}.out_chim_mapped.sam",
            merged_sam = DATA_DIR + "/mapped/noncoding_rna_chim/{sample}.merged.sam",
            out_full = DATA_DIR + "/mapped/noncoding_rna_single/{sample}.out_full_singlereads.sam",
            out_candidates = DATA_DIR + "/mapped/noncoding_rna_candidates/{sample}_out_chim_softclip.fastq.gz",
            smallrna_whole_bam = DATA_DIR + "/mapped/noncoding_rna_chim/{sample}.out_chim.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/disambiguate_and_filter.log"
    threads: 20
    resources:  mem = 50
    params: read_len = 100, 
            map_perc = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_perc_small"].min(),
            map_score = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_score_small"].min(),
    conda:  "wraps/disambiguate_and_filter/env.yaml"
    script: "wraps/disambiguate_and_filter/script.py"

rule alignment_softclip:
    input:  fastq = DATA_DIR + "/mapped/noncoding_rna_candidates/{sample}_out_chim_softclip.fastq.gz",
            genome = DATA_DIR + REF,
            index = DATA_DIR + "/index/STAR/SAindex",
    output: bam = DATA_DIR + "/mapped/softclip_genomic_noncoding/{sample}.hyb_candidates.sam",
    log:    run = DATA_DIR + "/sample_logs/{sample}/alignment_softclip.log"
    threads: 20
    resources:  mem = 34
    params: prefix = DATA_DIR + "/mapped/softclip_genomic_noncoding/{sample}.",
            map_perc = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_perc_softclip"].min(),
            map_score = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_score_softclip"].min(),
    conda:  "wraps/alignment_softclip/env.yaml"
    script: "wraps/alignment_softclip/script.py"

rule modify_bedgraph:
    input:  bed = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm.bedgraph"
    output: mod_bed = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm_modified.bedgraph"
    log:    run = DATA_DIR + "/sample_logs/{sample}/modify_bedgraph.log"
    threads: 20
    resources:  mem = 34
    conda:  "wraps/modify_bedgraph/env.yaml"
    script: "wraps/modify_bedgraph/script.py"

rule get_pairs:
    input:  smallrna_bam = DATA_DIR + "/mapped/noncoding_rna_chim/{sample}.out_chim_mapped.sam",
            rna_bam = DATA_DIR + "/mapped/softclip_genomic_noncoding/{sample}.hyb_candidates.sam",
            mod_rna_bed = DATA_DIR + "/mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm_modified.bedgraph",
            repeatmasker_bed = DATA_DIR + REPEATMASKER,
            smallrna_whole_bam = DATA_DIR + "/mapped/noncoding_rna_chim/{sample}.out_chim.sam",
            transcripts = DATA_DIR + TRANSCRIPTS,
            features = DATA_DIR + FEATURES,
            mirna_db = DATA_DIR + REF_MIR
    output: hyb_pairs = DATA_DIR + "/hyb_pairs/{sample}.hyb_pairs_deduplicated.tsv"
    log:    run = DATA_DIR + "/sample_logs/{sample}/get_pairs.log"
    threads: 1
    resources:  mem = 34
    params: prefix = DATA_DIR + "/hyb_pairs/cofold_pairs/",
            is_umi = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"is_umi"].min(),
    conda:  "wraps/get_pairs/env.yaml"
    script: "wraps/get_pairs/script.py"

rule filter_and_collapse:
    input:  hyb_pairs = DATA_DIR + "/hyb_pairs/{sample}.hyb_pairs_deduplicated.tsv",
            mirna_db = DATA_DIR + REF_MIR,
            trna_db = DATA_DIR + REF_TRNA,
            rrna_db = DATA_DIR + REF_RRNA,
            vaultrna_db = DATA_DIR + REF_VAULTRNA,
            snorna_db = DATA_DIR + REF_SNORNA,
            yrna_db = DATA_DIR + REF_YRNA,
            repeatmasker_bed = DATA_DIR + REPEATMASKER,
            transcripts = DATA_DIR + TRANSCRIPTS,
            features = DATA_DIR + FEATURES,
    output: pair_collapsed_unique = DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_unique.tsv",
            pair_collapsed_ambiguous = DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_ambiguous.tsv",          
    log:    run = DATA_DIR + "/sample_logs/{sample}/filter_and_collapse.log"
    threads: 1
    resources:  mem = 34
    params: prefix = DATA_DIR + "/hyb_pairs/cofold_pairs_collapsed/",
            is_umi = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"is_umi"].min(),
    conda:  "wraps/filter_and_collapse/env.yaml"
    script: "wraps/filter_and_collapse/script.py"


rule merge_replicates:
    input:  hyb_pairs = expand(DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_{{type}}.tsv", sample = SAMPLES),
            repeatmasker_bed = DATA_DIR + REPEATMASKER,
            transcripts = DATA_DIR + TRANSCRIPTS,
            features = DATA_DIR + FEATURES,
    output: merged = DATA_DIR + "/hyb_pairs/Merged.hybrids_deduplicated_filtered_collapsed_{type}.tsv",        
    log:    run = DATA_DIR + "/sample_logs/Merged_{type}/merge_replicates.log"
    threads: 1
    resources:  mem = 34
    params: prefix = DATA_DIR + "/hyb_pairs/cofold_pairs_collapsed_merged_{type}/all_types/",
            is_umi = cfg["is_umi"][1],
            is_ambiguous =  "{type}" 
    conda:  "wraps/merge_replicates/env.yaml"
    script: "wraps/merge_replicates/script.py"

rule unify_length:
    input:  hyb_pairs = DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_{type}.tsv",
            repeatmasker_bed = DATA_DIR + REPEATMASKER,
            transcripts = DATA_DIR + TRANSCRIPTS,
            features = DATA_DIR + FEATURES,
    output: pair_collapsed_unified = DATA_DIR + "/hyb_pairs/{sample}.hybrids_deduplicated_filtered_collapsed_unified_length_all_types_{type}_high_confidence.tsv",          
    log:    run = DATA_DIR + "/sample_logs/{sample}_{type}/unify_length.log"
    threads: 1
    resources:  mem = 34
    params: prefix = DATA_DIR + "/hyb_pairs/cofold_pairs_collapsed_unified_{type}/",
            is_umi = cfg["is_umi"][1],
            merged = "TRUE" if len(SAMPLES) > 1 else "FALSE",
            is_ambiguous =  "{type}"
    conda:  "wraps/unify_length/env.yaml"
    script: "wraps/unify_length/script.py"
