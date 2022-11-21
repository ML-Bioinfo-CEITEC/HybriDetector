######################################
# wrapper for rule: filter_single_genomic
######################################
import os
import sys
import math
import subprocess
import glob
from snakemake.shell import shell

UMITOOLS = "umi_tols"
SAMTOOLS = "samtools"
BAMCOVERAGE = "bedGraphToBigWig"

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: filter_single_genomic \n##\n")
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

command = "mkdir -p "+os.path.dirname(snakemake.params.prefix_non_genomic)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#keep header of alignment file
command = SAMTOOLS + " view -H -@ " + str(snakemake.threads) +" "+ snakemake.input.bam + " > " + snakemake.input.bam.replace(".to_genome.bam",".unique.sam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# get only genomic reads, mapped, qmap == 255 (unique), ratio alignment length/ matched part > 0.85 and no overlap junction (N) - 85% gives us maximum of 10bp softclip, which even though, it will be chimera, it cannot be determined at all
command = SAMTOOLS + " view -@ " + str(snakemake.threads) +" "+ snakemake.input.bam + " | awk ' {{ patsplit($6, arr, /([0-9]+M)/); total=0; for(i=1; i <= length(arr); ++i) {{total = total+arr[i]}}; if(($5==255) && ($6!~/N/) && total/length($10) > 0.85 ) print }}' - >> " + snakemake.input.bam.replace(".to_genome.bam",".unique.sam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#remove all alignments in mirna regions
command = SAMTOOLS + " view -b -h -@ " + str(snakemake.threads) +" "+ snakemake.input.bam.replace(".to_genome.bam",".unique.sam") + " -U " + snakemake.output.bam + " -L " + snakemake.input.bed + " > " + snakemake.input.bam.replace(".to_genome.bam",".unique_mirnaloc.bam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#get readnames of single genomic 
command = SAMTOOLS + " view -@ " + str(snakemake.threads) +" "+ snakemake.output.bam + " | cut -f1 > " + snakemake.output.bam.replace(".bam",".readids")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# here be carefull about grep, if there is . in read name, grep is extremely exhaustive and take a lot of memory unlesyou use -F
#second thing is that if we have order number in read names,it is highly probable to filter much more reads than we wanna, therefore use -w as "exact word", 
# but again be carefull, because the match muse be exact, nothing more in the begining or the end of word (apart space)
command = SAMTOOLS + " view -h -@ " + str(snakemake.threads) +" "+ snakemake.input.bam + " | grep -vwFf " +snakemake.output.bam.replace(".bam",".readids")+ " | " +SAMTOOLS + " view -bS -@ " + str(snakemake.threads) +" -o " + snakemake.output.bam_non_genomic + " -"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)    

command = SAMTOOLS + " bam2fq -@ " + str(snakemake.threads) +" "+ snakemake.output.bam_non_genomic + " > " +snakemake.output.bam_non_genomic.replace(".bam",".fastq")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "gzip "+snakemake.output.fastq.replace(".gz","")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS +" index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam 
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS +" index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam_non_genomic 
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


