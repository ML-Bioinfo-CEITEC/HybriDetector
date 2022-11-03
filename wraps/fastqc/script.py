######################################
# wrapper for rule: raw_fastq_qc
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: raw_fastq_qc \n##\n")
f.close()

TOOL = "fastqc"

shell.executable("/bin/bash")

version = str(subprocess.Popen(TOOL+" --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p `dirname "+snakemake.output.html+"`"+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = TOOL+" -o "+snakemake.params.prefix+" --threads "+str(snakemake.threads)+" "+snakemake.input.reads+" >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
