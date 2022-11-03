#############################################################
# wrapper for rule: modify_bedgraph
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: modify_bedgraph \n##\n")
f.close()

command = "$CONDA_PREFIX/bin/Rscript "+os.path.abspath(os.path.dirname(__file__))+"/modify_bed.R "+\
            snakemake.input.bed+ " " +snakemake.wildcards.sample + " >> " + snakemake.log.run + " 2>&1 "

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
