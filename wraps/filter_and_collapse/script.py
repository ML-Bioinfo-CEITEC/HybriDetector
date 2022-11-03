#############################################################
# wrapper for rule: filter_and_collapse
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: filter_and_collapse \n##\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix)+ "/" + snakemake.wildcards.sample + " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$CONDA_PREFIX/bin/Rscript " + os.path.abspath(os.path.dirname(__file__))+"/filter_and_collapse.R " + \
            snakemake.input.hyb_pairs + " " + snakemake.wildcards.sample + " " + snakemake.input.repeatmasker_bed + \
            " " + snakemake.input.transcripts + " " + snakemake.input.features +  " " + snakemake.input.mirna_db + \
            " " + snakemake.input.trna_db + " " + snakemake.input.rrna_db + " " + snakemake.input.vaultrna_db + \
            " " + snakemake.input.snorna_db + " " + snakemake.input.yrna_db + " " + snakemake.params.is_umi +" >> " + snakemake.log.run + " 2>&1 "

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
