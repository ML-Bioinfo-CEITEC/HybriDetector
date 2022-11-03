library(data.table)
run_all <- function(args){
	bedgraph <- args[1]
	sample <- args[2]
	rnabed <- fread(bedgraph)
	rnabed <- rnabed[V4 != 0]
	rnabed_vec <- numeric(nrow(rnabed) * 2)
	rnabed_vec[seq(1,length(rnabed_vec),by = 2)] <- rnabed$V2
	rnabed_vec[seq(2,length(rnabed_vec),by = 2)] <- rnabed$V3
	rnabed[,connect := diff(rnabed_vec)[seq(2,length(rnabed_vec),by = 2)] < 20]

	for(i in which(rnabed$connect == TRUE)){
	  set(rnabed,i+1L,"V2",rnabed$V2[i])
	  if(rnabed[i]$V4 > rnabed[i+1]$V4 ){
	    set(rnabed,i+1L,"V4",rnabed$V4[i])
	  }
	}

	#genomic mapping with any coverage
	names(rnabed) <- c("seqnames","start","end","cov","connect")
	rnabed <- rnabed[connect != TRUE]
	rnabed$connect <- NULL
	write.table(rnabed,paste0(c("mapped/genomic/",sample,".unique_genomic_dedup.without_norm.mapped_reads_norm_modified.bedgraph"),collapse = ""),col.names = F,row.names = F,quote = F,sep = "\t")

}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)