######################################
# wrapper for rule: alignment_single_genomic
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
f.write("\n##\n## RULE: alignment_SE_RNA \n##\n")
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

# to calculate map perc and map score automatically from read length - this will force minimum length of noncoding part to 16
perc = math.ceil(16/int(snakemake.params.read_length)*100)
command = STAR+" --runMode alignReads --runThreadN " + str(snakemake.threads) + \
               " --genomeDir " + snakemake.input.index.replace("/SAindex","") + \
               " --readFilesIn " + str(snakemake.input.fastq)  + \
               " --readFilesCommand zcat" + \
               " --outFileNamePrefix " + snakemake.params.prefix + \
               " --outFilterMultimapNmax 20 --outSAMmultNmax 3 --limitBAMsortRAM " + str(snakemake.resources.mem)+"000000000" \
               " --outFilterMismatchNmax 999" + \
               " --outFilterMismatchNoverReadLmax 0.1" + \
               " --outFilterMismatchNoverLmax 0.1" + \
               " --alignMatesGapMax 1000000" + \
               " --outFilterMatchNmin 0 --outFilterScoreMinOverLread " + str(perc)+" --outFilterMatchNminOverLread " + str(perc)+ \
               " --outSAMheaderHD @HD VN:1.4 SO:coordinate" + \
               " --outSAMattrRGline ID:"+str(snakemake.wildcards.sample)+" PL:Illumina PU:"+str(snakemake.wildcards.sample)+" SM:"+str(snakemake.wildcards.sample) + \
               " --outSAMunmapped Within --outSAMattributes All " + \
               " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS + " view -H -@ "+str(snakemake.threads)+ " " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam > " + snakemake.output.bam.replace(".mapped.filtered.sam","_header.sam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS + " view -F 276 -@ "+str(snakemake.threads)+ " " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam -o " + snakemake.output.bam.replace(".filtered.sam",".sam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#filter out reads shorter than 17 bp
command = "awk ' {{ patsplit($6, arr, /([0-9]+M)/); total=0; for(i=1; i <= length(arr); ++i) {{total = total+arr[i]}}; if((length($10) > 17) && ($6!~/N/) && total > 15 ) print }}' "+ snakemake.output.bam.replace(".filtered.sam",".sam") + " > " + snakemake.output.bam
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#convert to bam
command = "cat " + snakemake.output.bam.replace(".filtered.sam",".sam") + " >> " +snakemake.output.bam.replace(".mapped.filtered.sam","_header.sam") 
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS + " view -bh -@ "+str(snakemake.threads)+ " " + snakemake.output.bam.replace(".mapped.filtered.sam","_header.sam") +" > " + snakemake.output.bam.replace(".sam",".bam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS +" index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam.replace(".sam",".bam")
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

