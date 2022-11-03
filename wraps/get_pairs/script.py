#############################################################
# wrapper for rule: get_pairs
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: get_pairs \n##\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix)+ "/" + snakemake.wildcards.sample + " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$CONDA_PREFIX/bin/Rscript " + os.path.abspath(os.path.dirname(__file__))+"/get_pairs.R " + \
            snakemake.input.smallrna_bam + " " + snakemake.wildcards.sample + " " + snakemake.input.rna_bam + \
            " " + snakemake.input.mod_rna_bed + " " + snakemake.input.repeatmasker_bed +  " " + snakemake.input.smallrna_whole_bam + \
            " " + snakemake.input.transcripts + " " + snakemake.input.features + " " + snakemake.input.mirna_db + " " + snakemake.params.is_umi + " >> " + snakemake.log.run + " 2>&1 "

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
