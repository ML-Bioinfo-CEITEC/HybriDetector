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

REF_MIR = "DBs/mirna/hsa_mature.fa"
REF_TRNA = "DBs/tRNA_final.fa"
REF_RRNA = "DBs/rrna/rRNA_ensembl_dedup_final.fa"
REF_SNORNA = "DBs/snorna/snorna_dashr_dedup_final.fa"
REF_TRNA_ENS = "DBs/trna/ensembl/tRNA_ensembl_dedup_final.fa"
REF_TRNA_GTRNADB = "DBs/trna/gtrnadb/tRNA_gtrnadb_dedup_final.fa"
REF_TRNA_NCBI = "DBs/trna/ncbi/tRNA_ncbi_dedup_final.fa"
REF_TRNA_UCSC= "DBs/trna/ucsc/tRNA_ucsc_dedup_final.fa"
REF_VAULTRNA = "DBs/vaultrna/vaultRNA_dedup_final.fa"
REF_YRNA = "DBs/YRNA_final.fa"
REF_YRNA_ENS = "DBs/yrna/ensembl/YRNA_ensembl_dedup_final.fa"
REF_YRNA_UCSC = "DBs/yrna/ucsc/YRNA_ucsc_dedup_final.fa"
NON_CODING_RNA_BED = "DBs/non_coding_rna_sorted_merged.bed"
REPEATMASKER = "DBs/repeatmasker.bed"
TRANSCRIPTS = "DBs/GRCh38-p10_transcripts_of_longest_transcripts.gtf"
FEATURES = "DBs/GRCh38-p10_longest_transcripts_subtracted_final.gtf" 


#DATA_DIR = "/mnt/ssd/ssd_3/temp/vasek/CLASH"
SAMPLES = cfg["Sample"]

wildcard_constraints:
    type="(unique|ambiguous)"



def all_inputs(wildcards):
    if len(SAMPLES) > 1:
        return expand("hyb_pairs/{sample}.unified_length_all_types_{type}_high_confidence.tsv", sample = "Merged", type = ["unique","ambiguous"])
    else: 
        return expand("hyb_pairs/{sample}.unified_length_all_types_{type}_high_confidence.tsv", sample = SAMPLES, type = ["unique","ambiguous"])

rule all:
    input:
        pair_collapsed_unified = all_inputs

rule prepare_ref:
    input: 
        features = FEATURES + ".gz",
        repeatmasker = REPEATMASKER + ".gz"
    output:
        genome = "DBs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        features = "DBs/GRCh38-p10_longest_transcripts_subtracted_final.gtf",
        repeatmasker = "DBs/repeatmasker.bed" 
    log:
        run = "sample_logs/prepare_ref.log"
    params:
        prefix = "DBs/"
    threads:    1
    conda:  "wraps/prepare_ref/env.yaml"
    script: "wraps/prepare_ref/script.py"

REF = "DBs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

rule fastqc:
	input:	
		reads = "data/{sample}.fastq.gz"
	output:
		html = "raw_fastq_qc/{sample}.fastqc.html"
	log:
		run = "sample_logs/{sample}/raw_fastq_qc.log"
	params:	
		prefix = "raw_fastq_qc/",
		html = "raw_fastq_qc/{sample}_fastqc.html"
	threads:	1
	conda:	"wraps/fastqc/env.yaml"
	script:	"wraps/fastqc/script.py"


rule STAR_gen_index:
    input:  gen = REF,       
    output: SAindex = "index/STAR/SAindex",
    params: dir = "index/STAR",
            log = "index/STAR/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR.indexation_run.log",
    threads:    40
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_mirna_gen_index:
    input:  gen = REF_MIR,       
    output: SAindex = "index/STAR_mirna/SAindex",
    params: dir = "index/STAR_mirna",
            log = "index/STAR_mirna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_mirna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_rrna_gen_index:
    input:  gen = REF_RRNA,       
    output: SAindex = "index/STAR_rrna/SAindex",
    params: dir = "index/STAR_rrna",
            log = "index/STAR_rrna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_rrna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_snorna_gen_index:
    input:  gen = REF_SNORNA,       
    output: SAindex = "index/STAR_snorna/SAindex",
    params: dir = "index/STAR_snorna",
            log = "index/STAR_snorna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_snorna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ens_gen_index:
    input:  gen = REF_TRNA_ENS,       
    output: SAindex = "index/STAR_trna_ens/SAindex",
    params: dir = "index/STAR_trna_ens",
            log = "index/STAR_trna_ens/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ens.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_gtrnadb_gen_index:
    input:  gen = REF_TRNA_GTRNADB,       
    output: SAindex = "index/STAR_trna_gtrnadb/SAindex",
    params: dir = "index/STAR_trna_gtrnadb",
            log = "index/STAR_trna_gtrnadb/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_gtrnadb.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ncbi_gen_index:
    input:  gen = REF_TRNA_NCBI,       
    output: SAindex = "index/STAR_trna_ncbi/SAindex",
    params: dir = "index/STAR_trna_ncbi",
            log = "index/STAR_trna_ncbi/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ncbi.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_trna_ucsc_gen_index:
    input:  gen = REF_TRNA_UCSC,       
    output: SAindex = "index/STAR_trna_ucsc/SAindex",
    params: dir = "index/STAR_trna_ucsc",
            log = "index/STAR_trna_ucsc/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_trna_ucsc.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_vaultrna_gen_index:
    input:  gen = REF_VAULTRNA,       
    output: SAindex = "index/STAR_vaultrna/SAindex",
    params: dir = "index/STAR_vaultrna",
            log = "index/STAR_vaultrna/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_vaultrna.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_yrna_ens_index:
    input:  gen = REF_YRNA_ENS,       
    output: SAindex = "index/STAR_yrna_ens/SAindex",
    params: dir = "index/STAR_yrna_ens",
            log = "index/STAR_yrna_ens/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_yrna_ens.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule STAR_yrna_ucsc_gen_index:
    input:  gen = REF_YRNA_UCSC,       
    output: SAindex = "index/STAR_yrna_ucsc/SAindex",
    params: dir = "index/STAR_yrna_ucsc",
            log = "index/STAR_yrna_ucsc/Log.out",
    resources:  mem = 200
    log:    run = "index/STAR_yrna_ucsc.indexation_run.log",
    threads:    30
    conda:  "wraps/STAR_gen_index/env.yaml"
    script: "wraps/STAR_gen_index/script.py"

rule alignment_single_genomic:
    input:  fastq = "data/{sample}.fastq.gz",
            genome = REF,
            index = "index/STAR/SAindex",
    output: bam = "mapped/to_genome/{sample}.to_genome.bam",
            bai = "mapped/to_genome/{sample}.to_genome.bam.bai",
    log:    run = "sample_logs/{sample}/alignment.log"
    threads: 20
    resources:  mem = 34
    params: prefix = "mapped/to_genome/{sample}/{sample}.",
            read_len = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min(), # Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html
            map_perc = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_perc_single_genomic"].min(),
            snakemake_dir = workflow.basedir + "/../"
    conda:  "wraps/alignment/env.yaml"
    script: "wraps/alignment/script.py"

rule filter_single_genomic:
    input:  bam = "mapped/to_genome/{sample}.to_genome.bam",
            bai = "mapped/to_genome/{sample}.to_genome.bam.bai",
            bed = NON_CODING_RNA_BED
    output: bam = "mapped/genomic/{sample}.unique_genomic.bam",
            bam_non_genomic = "mapped/non_genomic/{sample}.unique.removedGenomic.bam",
            fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz"
    log:    run = "sample_logs/{sample}/filter_single_genomic.log"
    threads: 20
    resources:  mem = 34
    params: prefix = "mapped/genomic/{sample}.",
            prefix_non_genomic = "mapped/non_genomic/{sample}."          
    conda:  "wraps/filter_single_genomic/env.yaml"
    script: "wraps/filter_single_genomic/script.py"

rule deduplicate_single_genomic:
    input: bam = "mapped/genomic/{sample}.unique_genomic.bam"
    output: bam = "mapped/genomic/{sample}.unique_genomic_dedup.bam",
            bed = "mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm.bedgraph"
    log:    run = "sample_logs/{sample}/deduplicate_single_genomic.log"
    threads: 20
    resources:  mem = 34
    params:
        is_umi = cfg["is_umi"][0],
    conda:  "wraps/deduplicate_single_genomic/env.yaml"
    script: "wraps/deduplicate_single_genomic/script.py"

rule alignment_mirna:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_MIR,
            index = "index/STAR_mirna/SAindex",
    output: bam = "mapped/non_genomic/mirna_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_mirna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/mirna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_rrna:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_RRNA,
            index = "index/STAR_rrna/SAindex",
    output: bam = "mapped/non_genomic/rrna_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_rrna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/rrna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_snorna:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_SNORNA,
            index = "index/STAR_snorna/SAindex",
    output: bam = "mapped/non_genomic/snorna_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_snorna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/snorna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ens:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_TRNA_ENS,
            index = "index/STAR_trna_ens/SAindex",
    output: bam = "mapped/non_genomic/trna_ens_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_trna_ens.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/trna_ens_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_gtrnadb:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_TRNA_GTRNADB,
            index = "index/STAR_trna_gtrnadb/SAindex",
    output: bam = "mapped/non_genomic/trna_gtrnadb_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_trna_gtrnadb.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/trna_gtrnadb_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ncbi:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_TRNA_NCBI,
            index = "index/STAR_trna_ncbi/SAindex",
    output: bam = "mapped/non_genomic/trna_ncbi_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_trna_ncbi.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/trna_ncbi_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_trna_ucsc:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_TRNA_UCSC,
            index = "index/STAR_trna_ucsc/SAindex",
    output: bam = "mapped/non_genomic/trna_ucsc_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_trna_ucsc.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/trna_ucsc_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_vaultrna:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_VAULTRNA,
            index = "index/STAR_vaultrna/SAindex",
    output: bam = "mapped/non_genomic/vaultrna_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_vaultrna.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/vaultrna_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_yrna_ens:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_YRNA_ENS,
            index = "index/STAR_yrna_ens/SAindex",
    output: bam = "mapped/non_genomic/yrna_ens_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_yrna_ens.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/yrna_ens_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule alignment_yrna_ucsc:
    input:  fastq = "mapped/non_genomic/{sample}.unique.removedGenomic.fastq.gz",
            genome = REF_YRNA_UCSC,
            index = "index/STAR_yrna_ucsc/SAindex",
    output: bam = "mapped/non_genomic/yrna_ucsc_alignment/{sample}.mapped.filtered.sam",
    log:    run = "sample_logs/{sample}/alignment_yrna_ucsc.log"
    threads: 20
    resources:  mem = 50
    params: prefix = "mapped/non_genomic/yrna_ucsc_alignment/{sample}.",
            read_len = 100, 
            read_length = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"read_length"].min()
    conda:  "wraps/alignment_small/env.yaml"
    script: "wraps/alignment_small/script.py"

rule disambiguate_and_filter:
    input:  mirna_sam = "mapped/non_genomic/mirna_alignment/{sample}.mapped.filtered.sam",
            rrna_sam = "mapped/non_genomic/rrna_alignment/{sample}.mapped.filtered.sam",
            snorna_sam = "mapped/non_genomic/snorna_alignment/{sample}.mapped.filtered.sam",
            trna_ens_sam = "mapped/non_genomic/trna_ens_alignment/{sample}.mapped.filtered.sam",
            trna_gtrnadb_sam = "mapped/non_genomic/trna_gtrnadb_alignment/{sample}.mapped.filtered.sam",
            trna_ncbi_sam = "mapped/non_genomic/trna_ncbi_alignment/{sample}.mapped.filtered.sam",
            trna_ucsc_sam = "mapped/non_genomic/trna_ucsc_alignment/{sample}.mapped.filtered.sam",
            vaultrna_sam = "mapped/non_genomic/vaultrna_alignment/{sample}.mapped.filtered.sam",
            yrna_ens_sam = "mapped/non_genomic/yrna_ens_alignment/{sample}.mapped.filtered.sam",
            yrna_ucsc_sam = "mapped/non_genomic/yrna_ucsc_alignment/{sample}.mapped.filtered.sam",
    output: bam = "mapped/noncoding_rna_chim/{sample}.out_chim_mapped.sam",
            merged_sam = "mapped/noncoding_rna_chim/{sample}.merged.sam",
            out_full = "mapped/noncoding_rna_single/{sample}.out_full_singlereads.sam",
            out_candidates = "mapped/noncoding_rna_candidates/{sample}_out_chim_softclip.fastq.gz",
            smallrna_whole_bam = "mapped/noncoding_rna_chim/{sample}.out_chim.sam",
    log:    run = "sample_logs/{sample}/disambiguate_and_filter.log"
    threads: 20
    resources:  mem = 50
    conda:  "wraps/disambiguate_and_filter/env.yaml"
    script: "wraps/disambiguate_and_filter/script.py"

rule alignment_softclip:
    input:  fastq = "mapped/noncoding_rna_candidates/{sample}_out_chim_softclip.fastq.gz",
            genome = REF,
            index = "index/STAR/SAindex",
    output: bam = "mapped/softclip_genomic_noncoding/{sample}.hyb_candidates.sam",
    log:    run = "sample_logs/{sample}/alignment_softclip.log"
    threads: 20
    resources:  mem = 34
    params: prefix = "mapped/softclip_genomic_noncoding/{sample}.",
            map_perc = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"map_perc_softclip"].min(),
    conda:  "wraps/alignment_softclip/env.yaml"
    script: "wraps/alignment_softclip/script.py"

rule modify_bedgraph:
    input:  bed = "mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm.bedgraph"
    output: mod_bed = "mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm_modified.bedgraph"
    log:    run = "sample_logs/{sample}/modify_bedgraph.log"
    threads: 20
    resources:  mem = 34
    conda:  "wraps/modify_bedgraph/env.yaml"
    script: "wraps/modify_bedgraph/script.py"

rule get_pairs:
    input:  smallrna_bam = "mapped/noncoding_rna_chim/{sample}.out_chim_mapped.sam",
            rna_bam = "mapped/softclip_genomic_noncoding/{sample}.hyb_candidates.sam",
            mod_rna_bed = "mapped/genomic/{sample}.unique_genomic_dedup.without_norm.mapped_reads_norm_modified.bedgraph",
            repeatmasker_bed = REPEATMASKER,
            smallrna_whole_bam = "mapped/noncoding_rna_chim/{sample}.out_chim.sam",
            transcripts = TRANSCRIPTS,
            features = FEATURES,
            mirna_db = REF_MIR
    output: hyb_pairs = "hyb_pairs/{sample}.hyb_pairs_deduplicated.tsv"
    log:    run = "sample_logs/{sample}/get_pairs.log"
    threads: 1
    resources:  mem = 34
    params: prefix = "hyb_pairs/cofold_pairs/",
            is_umi = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"is_umi"].min(),
    conda:  "wraps/get_pairs/env.yaml"
    script: "wraps/get_pairs/script.py"

rule filter_and_collapse:
    input:  hyb_pairs = "hyb_pairs/{sample}.hyb_pairs_deduplicated.tsv",
            mirna_db = REF_MIR,
            trna_db = REF_TRNA,
            rrna_db = REF_RRNA,
            vaultrna_db = REF_VAULTRNA,
            snorna_db = REF_SNORNA,
            yrna_db = REF_YRNA,
            repeatmasker_bed = REPEATMASKER,
            transcripts = TRANSCRIPTS,
            features = FEATURES,
    output: pair_collapsed_unique = "hyb_pairs/{sample}.unique.tsv",
            pair_collapsed_ambiguous = "hyb_pairs/{sample}.ambiguous.tsv",          
    log:    run = "sample_logs/{sample}/filter_and_collapse.log"
    threads: 1
    resources:  mem = 34
    params: prefix = "hyb_pairs/cofold_pairs_collapsed/",
            is_umi = lambda wildcards: cfg.loc[SAMPLES == wildcards.sample,"is_umi"].min(),
    conda:  "wraps/filter_and_collapse/env.yaml"
    script: "wraps/filter_and_collapse/script.py"


rule merge_replicates:
    input:  hyb_pairs = expand("hyb_pairs/{sample}.{{type}}.tsv", sample = SAMPLES),
            repeatmasker_bed = REPEATMASKER,
            transcripts = TRANSCRIPTS,
            features = FEATURES,
    output: merged = "hyb_pairs/Merged.{type}.tsv",        
    log:    run = "sample_logs/Merged_{type}/merge_replicates.log"
    threads: 1
    resources:  mem = 34
    params: prefix = "hyb_pairs/cofold_pairs_collapsed_merged_{type}/all_types/",
            is_umi = cfg["is_umi"][0],
            is_ambiguous =  "{type}" 
    conda:  "wraps/merge_replicates/env.yaml"
    script: "wraps/merge_replicates/script.py"

rule unify_length:
    input:  hyb_pairs = "hyb_pairs/{sample}.{type}.tsv",
            repeatmasker_bed = REPEATMASKER,
            transcripts = TRANSCRIPTS,
            features = FEATURES,
    output: pair_collapsed_unified = "hyb_pairs/{sample}._unified_length_all_types_{type}_high_confidence.tsv",          
    log:    run = "sample_logs/{sample}_{type}/unify_length.log"
    threads: 1
    resources:  mem = 34
    params: prefix = "hyb_pairs/cofold_pairs_collapsed_unified_{type}/",
            is_umi = cfg["is_umi"][0],
            merged = "TRUE" if len(SAMPLES) > 1 else "FALSE",
            is_ambiguous =  "{type}"
    conda:  "wraps/unify_length/env.yaml"
    script: "wraps/unify_length/script.py"
