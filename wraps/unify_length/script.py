#############################################################
# wrapper for rule: unify_length
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: unify_length \n##\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix)+ "/" + snakemake.wildcards.sample + " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

print(snakemake.input.hyb_pairs)
command = "$CONDA_PREFIX/bin/Rscript " + os.path.abspath(os.path.dirname(__file__))+"/unify_length.R " + \
            snakemake.input.hyb_pairs + " " + snakemake.wildcards.sample + " " + snakemake.input.repeatmasker_bed + \
            " " + snakemake.input.transcripts + " " + snakemake.input.features + " " + snakemake.params.is_umi + " " + \
            snakemake.params.merged+ " " + snakemake.params.is_ambiguous + " >> " + snakemake.log.run + " 2>&1 "

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
