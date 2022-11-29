######################################
# wrapper for rule: alignment_softclip
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

STAR = "STAR"
SAMTOOLS = "samtools"

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: alignment_softclip \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = STAR+" --runMode alignReads --runThreadN " + str(snakemake.threads) + \
               " --genomeDir " + snakemake.input.index.replace("/SAindex","") + \
               " --readFilesIn " + str(snakemake.input.fastq)  + \
               " --readFilesCommand zcat" + \
               " --outFileNamePrefix " + snakemake.params.prefix + \
               " --outFilterMultimapNmax 20" + \
               " --outFilterMismatchNmax 999" + \
               " --outFilterMismatchNoverReadLmax 0.1" + \
               " --outFilterMismatchNoverLmax 0.1" + \
               " --alignMatesGapMax 1000000" + \
               " --outFilterMatchNmin 0 --outFilterScoreMinOverLread " + str(snakemake.params.map_perc)+" --outFilterMatchNminOverLread " + str(snakemake.params.map_perc)+ \
               " --outSAMheaderHD @HD VN:1.4 SO:coordinate" + \
               " --outSAMattrRGline ID:"+str(snakemake.wildcards.sample)+" PL:Illumina PU:"+str(snakemake.wildcards.sample)+" SM:"+str(snakemake.wildcards.sample) + \
               " --outSAMunmapped Within --outSAMattributes All " + \
               " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS + " view -F 4 -@ " + str(snakemake.threads) + " " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam | awk '{{ if(($5==255) && ($2!=\"*\") && ($6!~/N/) ) print }}' - " + " > " + snakemake.output.bam
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
