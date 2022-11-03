######################################
# wrapper for rule: disambiguate_and_filter
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

STAR = "STAR"
SAMTOOLS = "samtools"
BEDGRAPH2BIGWIG = "bedGraphToBigWig"

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: disambiguate_and_filter \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "cat " + snakemake.input.mirna_sam + " " + snakemake.input.rrna_sam + " " + snakemake.input.snorna_sam + " " + snakemake.input.trna_ens_sam + " " + snakemake.input.trna_gtrnadb_sam + " " + snakemake.input.trna_ncbi_sam + " " + snakemake.input.trna_ucsc_sam + " " + snakemake.input.vaultrna_sam + " " + snakemake.input.yrna_ens_sam + " " + snakemake.input.yrna_ucsc_sam + " > " + snakemake.output.merged_sam 
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "sort -k1 --temporary-directory=/mnt/ssd/ssd_3/temp/vasek -S 40% --parallel=50 " + snakemake.output.merged_sam + " > " + snakemake.output.merged_sam.replace("merged.sam","merged_sorted.sam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "perl " + os.path.abspath(os.path.dirname(__file__))+"/sam_disambiguate.pl "+\
          " --verbose --ifile " + snakemake.output.merged_sam.replace("merged.sam","merged_sorted.sam")+\
          " --ofile_full " + snakemake.output.out_full+ \
          " --ofile_chim " + snakemake.output.smallrna_whole_bam + \
          " --ofile_chim_softclip " + snakemake.output.out_candidates.replace(".gz","") + \
          " --ofile_chim_mapped " + snakemake.output.merged_sam.replace("merged.sam","out_chim_mapped.sam")+ \
          " --ofile_chim_refref " + snakemake.output.merged_sam.replace("merged.sam","out_chim_refref.sam") +">> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "gzip " + snakemake.output.out_candidates.replace(".gz","")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)