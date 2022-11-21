######################################
# wrapper for rule: deduplicate_single_genomics
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

UMITOOLS = "umi_tools"
SAMTOOLS = "samtools"
BAMCOVERAGE = "bamCoverage"

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: deduplicate_single_genomics \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.is_umi == "FALSE":
    TOOL = "picard MarkDuplicates"
    SAMTOOLS = "samtools"

    version = str(subprocess.Popen(TOOL+" --version 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f = open(snakemake.log.run, 'at')
    f.write("## VERSION: Picard "+version+"\n")
    f.close()

    command = "export LD_BIND_NOW=1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = TOOL+" INPUT="+snakemake.input.bam+" OUTPUT="+snakemake.output.bam+" METRICS_FILE="+snakemake.output.bam.replace(".bam",".log")+" REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx"+str(snakemake.resources.mem)+"g 2>> "+snakemake.log.run+" "
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

else:
	command = UMITOOLS + " dedup -I " + snakemake.input.bam + " -S " + snakemake.output.bam + " --log " + snakemake.output.bam.replace(".bam",".log") + " --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 --spliced-is-unique --multimapping-detection-method=NH"
	f = open(snakemake.log.run, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)

command = SAMTOOLS +" index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam 
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = BAMCOVERAGE +" --normalizeUsing None -of bedgraph --binSize 1 -p "+str(snakemake.threads)+ " -b "+ snakemake.output.bam + " -o " + snakemake.output.bam.replace(".bam",".without_norm.bedgraph") + " >> " + snakemake.log.run + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

DIV= str(subprocess.Popen(SAMTOOLS+ " view -F 4 -@ " + str(snakemake.threads) +" "+ snakemake.output.bam +" | cut -f1 | sort -u -S20G --parallel=20 -T /mnt/ssd/ssd_1/tmp/ | wc -l | tr -d '\n'", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## div: "+DIV+"\n")
f.close()

command = "cat "+ snakemake.output.bam.replace(".bam",".without_norm.bedgraph") + " | awk -v divisor="+DIV+" 'BEGIN {{OFS=\"\\t\"}};{{$4=$4/(divisor/1000000); print $0;}}' > "+ snakemake.output.bed
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
