######################################
# wrapper for rule: STAR_gen_index
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

TOOL = "STAR"

help = subprocess.Popen("grep -v '>' " + snakemake.input.gen + " | wc -m",shell=True,stdout=subprocess.PIPE).communicate()[0]
STAR_GENOME_BASES_LOG = min(14,math.floor(math.log(float(int(help)),2)/2-1))

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'wt')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+snakemake.params.dir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

IS_COMPRESSED = str(subprocess.Popen("file " + snakemake.input.gen + " | grep 'compressed'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')

if IS_COMPRESSED:
	command = "gunzip -fk "+snakemake.input.gen+" >> "+snakemake.log.run+" 2>&1"
	f = open(snakemake.log.run, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)

command = TOOL+" --runMode genomeGenerate --runThreadN "+str(snakemake.threads)+" --limitGenomeGenerateRAM " + str(snakemake.resources.mem * 1000000000) +" --genomeDir "+snakemake.params.dir+" --genomeFastaFiles "+snakemake.input.gen.replace(".gz","")+ " --genomeSAindexNbases "+str(STAR_GENOME_BASES_LOG)+" >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -rf _STARtmp " + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
