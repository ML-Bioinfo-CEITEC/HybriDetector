library(data.table)
library(seqinr)
library(miRBaseConverter)
library(Biostrings)
library(stringi)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

run_all <- function(args){
	hybrids <- args[1]
	sample <- args[2]
    repeatmasker_bed <-args[3]
    transcripts <- args[4]
    features <- args[5]
    mirna_db <- args[6]
    trna_db <- args[7]
    rrna_db <- args[8]
    vaultrna_db <- args[9]
    snorna_db <- args[10]
    yrna_db <- args[11]
    is_umi <- args[12]

	pair_occurence_sequence <- fread(hybrids)

	##########################################
	#repeatmasker to remove 
	remove_repeatmasker <- c("Helitron1Na_Mam","(G)n")
	pair_occurence_sequence[,driver_length:= nchar(seq.m)]

	pair_occurence_sequence <- pair_occurence_sequence[order(chr.g,smallrna,pos.g)]
	 ###############################
	###collapse reads
	pair_occurence_sequence_filt <- pair_occurence_sequence[mrna_length == TRUE & mrna_aln == TRUE & !repeatmasker %in% remove_repeatmasker & !(A_perc > 50 | C_perc > 50 | G_perc > 50 | T_perc > 50 )& driver_length <= 30]
	pair_occurence_sequence_filt_print <- pair_occurence_sequence[mrna_length == TRUE & mrna_aln == TRUE & !repeatmasker %in% remove_repeatmasker & !(A_perc > 50 | C_perc > 50 | G_perc > 50 | T_perc > 50 & driver_length <= 30)]
	pair_occurence_sequence_filt$driver_length <- NULL
	#collapse all targets starting on the same position in genome, get the longest alignment across these targets and sum Ndups and Nunique
	# to deduplicate same targets with different UMIs
	pair_occurence_sequence_filt <- pair_occurence_sequence_filt[,.(seq.g.len=max(nchar(seq.g)),seq.m=max(seq.m),Ndups=sum(Ndups),Nunique=sum(Nunique)), by=.(chr.g,pos.g,flag.g,smallrna,smallRNA_fam,strand,forward_dir,smallrna_nomm,mrna_length,
                                                                                                                                                            cov,Nreads,mrna_aln,mrna_aln2,mrna_aln5,mrna_aln10,mrna_aln20)]


	####################################################
	#add non-coding RNA type
	pair_occurence_sequence_filt[, noncodingRNA_type := sapply(seq_along(chr.g),function(x) tail(unlist(tstrsplit(smallrna[x],"_")),n=1))]
	pair_occurence_sequence_filt[noncodingRNA_type == 10]$noncodingRNA_type <- "YRNA"
	pair_occurence_sequence_filt[noncodingRNA_type == "RNA"]$noncodingRNA_type <- "YRNA"
	pair_occurence_sequence_filt[grepl("hsa",noncodingRNA_type)]$noncodingRNA_type <- "miRNA"

	#add sequence of full noncoding RNA
	#####################################
	#mirna
	mirnadb <- read.fasta(mirna_db,as.string = TRUE)
	seq_name <- names(mirnadb)
	sequence <- paste(mirnadb)
	mirnaDB <- data.table(seq_name, sequence)
	mirnaDB[, sequence := toupper(sequence)]
	#trna
	trnadb <- read.fasta(trna_db,as.string = TRUE)
	seq_name <- names(trnadb)
	sequence <- paste(trnadb)
	trnadb <- data.table(seq_name, sequence)
	trnadb[, sequence := toupper(sequence)]
	#rrna
	rrnadb <- read.fasta(rrna_db,as.string = TRUE)
	seq_name <- names(rrnadb)
	sequence <- paste(rrnadb)
	rrnadb <- data.table(seq_name, sequence)
	rrnadb[, sequence := toupper(sequence)]
	#vaultrna
	vaultrnadb <- read.fasta(vaultrna_db,as.string = TRUE)
	seq_name <- names(vaultrnadb)
	sequence <- paste(vaultrnadb)
	vaultrnadb <- data.table(seq_name, sequence)
	vaultrnadb[, sequence := toupper(sequence)]
	#snorna
	snornadb <- read.fasta(snorna_db,as.string = TRUE)
	seq_name <- names(snornadb)
	sequence <- paste(snornadb)
	snornadb <- data.table(seq_name, sequence)
	snornadb[, sequence := toupper(sequence)]
	#Yrna
	yrnadb <- read.fasta(yrna_db,as.string = TRUE)
	seq_name <- names(yrnadb)
	sequence <- paste(yrnadb)
	yrnadb <- data.table(seq_name, sequence)
	yrnadb[, sequence := toupper(sequence)]

	wholedb <- rbind(mirnaDB,trnadb,rrnadb,vaultrnadb,snornadb,yrnadb)
	pair_occurence_sequence_filt <- merge(pair_occurence_sequence_filt,wholedb,by.x = "smallrna",by.y = "seq_name", all.x =T)
	names(pair_occurence_sequence_filt)[names(pair_occurence_sequence_filt) == 'sequence'] <- 'noncodingRNA_seq'

	#first collapse by position and mirna name - faster
	pair_occurence_sequence_filt <- pair_occurence_sequence_filt[order(chr.g,pos.g)]
	pair_occurence_sequence_filt[,new_start:=c(-21,pos.g[-length(pos.g)]), by=.(chr.g,smallrna)]
	pair_occurence_sequence_filt[,new_start:=fifelse(pos.g <= new_start+20, FALSE, TRUE)]
	pair_occurence_sequence_filt[,new_pos:=pos.g]
	for(i in which(! pair_occurence_sequence_filt$new_start)){
	if(i > 1)
	  set(pair_occurence_sequence_filt, i, "new_pos", pair_occurence_sequence_filt$new_pos[i-1L])
	}

	if (is_umi == TRUE){
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt[,.(seq.g.len = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),seq.g.len],
		                                                                       smallRNA_fam = smallRNA_fam[which.max(Nunique)],
		                                                                       forward_dir = as.logical(max(forward_dir)),
		                                                                       Ndups = sum(Ndups),
		                                                                       Nunique = sum(Nunique),
		                                                                       flag.g = flag.g[which.max(Nunique)],
		                                                                       strand = strand[which.max(Nunique)],
		                                                                       smallrna_nomm = smallrna_nomm[which.max(Nunique)], 
		                                                                       mrna_length = mrna_length[which.max(Nunique)],
		                                                                       cov = sum(cov),
		                                                                       Nreads = sum(Nreads),
		                                                                       mrna_aln = mrna_aln[which.max(Nunique)], 
		                                                                       mrna_aln2 = mrna_aln2[which.max(Nunique)], 
		                                                                       mrna_aln5 = mrna_aln5[which.max(Nunique)], 
		                                                                       mrna_aln10 = mrna_aln10[which.max(Nunique)],
		                                                                       mrna_aln20 = mrna_aln20[which.max(Nunique)],
		                                                                       seq.m = seq.m[which.max(Nunique)],
		                                                                       noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                       noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)]
		),
		by=.(chr.g,new_pos,smallrna)]
	} else {
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt[,.(seq.g.len = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),seq.g.len],
		                                                                       smallRNA_fam = smallRNA_fam[which.max(Ndups)],
		                                                                       forward_dir = as.logical(max(forward_dir)),
		                                                                       Ndups = sum(Ndups),
		                                                                       flag.g = flag.g[which.max(Ndups)],
		                                                                       strand = strand[which.max(Ndups)],
		                                                                       smallrna_nomm = smallrna_nomm[which.max(Ndups)], 
		                                                                       mrna_length = mrna_length[which.max(Ndups)],
		                                                                       cov = sum(cov),
		                                                                       Nreads = sum(Nreads),
		                                                                       mrna_aln = mrna_aln[which.max(Ndups)], 
		                                                                       mrna_aln2 = mrna_aln2[which.max(Ndups)], 
		                                                                       mrna_aln5 = mrna_aln5[which.max(Ndups)], 
		                                                                       mrna_aln10 = mrna_aln10[which.max(Ndups)],
		                                                                       mrna_aln20 = mrna_aln20[which.max(Ndups)],
		                                                                       seq.m = seq.m[which.max(Ndups)],
		                                                                       noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                       noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)]
		),
		by=.(chr.g,new_pos,smallrna)]
	}  
	pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[order(chr.g,new_pos)]

	#second collapse by position and smallrna sequence - faster
	if (is_umi == TRUE){
		pair_occurence_sequence_filt_collap_name <- pair_occurence_sequence_filt_collap[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                                                   seq.g.len = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),seq.g.len],
		                                                                                   smallRNA_fam = smallRNA_fam[which.max(Nunique)],
		                                                                                   forward_dir = as.logical(max(forward_dir)),
		                                                                                   Ndups = sum(Ndups),
		                                                                                   Nunique = sum(Nunique), 
		                                                                                   flag.g = flag.g[which.max(Nunique)],
		                                                                                   strand = strand[which.max(Nunique)],
		                                                                                   smallrna_nomm = smallrna_nomm[which.max(Nunique)], 
		                                                                                   mrna_length = mrna_length[which.max(Nunique)],
		                                                                                   cov = sum(cov),
		                                                                                   Nreads = sum(Nreads),
		                                                                                   mrna_aln = mrna_aln[which.max(Nunique)], 
		                                                                                   mrna_aln2 = mrna_aln2[which.max(Nunique)], 
		                                                                                   mrna_aln5 = mrna_aln5[which.max(Nunique)], 
		                                                                                   mrna_aln10 = mrna_aln10[which.max(Nunique)],
		                                                                                   mrna_aln20 = mrna_aln20[which.max(Nunique)],
		                                                                                   noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                                   noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)]
		),
		by=.(chr.g,new_pos,seq.m)]
	} else {
		pair_occurence_sequence_filt_collap_name <- pair_occurence_sequence_filt_collap[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                                                   seq.g.len = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),seq.g.len],
		                                                                                   smallRNA_fam = smallRNA_fam[which.max(Ndups)],
		                                                                                   forward_dir = as.logical(max(forward_dir)),
		                                                                                   Ndups = sum(Ndups),
		                                                                                   flag.g = flag.g[which.max(Ndups)],
		                                                                                   strand = strand[which.max(Ndups)],
		                                                                                   smallrna_nomm = smallrna_nomm[which.max(Ndups)], 
		                                                                                   mrna_length = mrna_length[which.max(Ndups)],
		                                                                                   cov = sum(cov),
		                                                                                   Nreads = sum(Nreads),
		                                                                                   mrna_aln = mrna_aln[which.max(Ndups)], 
		                                                                                   mrna_aln2 = mrna_aln2[which.max(Ndups)], 
		                                                                                   mrna_aln5 = mrna_aln5[which.max(Ndups)], 
		                                                                                   mrna_aln10 = mrna_aln10[which.max(Ndups)],
		                                                                                   mrna_aln20 = mrna_aln20[which.max(Ndups)],
		                                                                                   noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                                   noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)]
		),
		by=.(chr.g,new_pos,seq.m)]
	}  
	pair_occurence_sequence_filt_collap_name <- pair_occurence_sequence_filt_collap_name[order(chr.g,new_pos)]
	  
	clustering <- function(input_tab) {
		#perform sequence similarity clustering
		distance_noncoding <- adist(input_tab$seq.m,input_tab$seq.m)
		mean_len_mat_noncoding <- as.matrix(proxy::dist(nchar(input_tab$seq.m),method = function(x,y) (x+y)/2))
		diag(mean_len_mat_noncoding) <- nchar(input_tab$seq.m)
		input_tab[,noncoding_clust := cutree(hclust(as.dist(distance_noncoding/mean_len_mat_noncoding)),h = 0.3)]
		#collapse by identified clusters
		if (is_umi == TRUE){
		  pair_occurence_sequence_filt_collap_clust <- input_tab[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                            flag.g = flag.g[which.max(Nunique)],
		                                                            strand = strand[which.max(Nunique)],
		                                                            smallRNA_fam = smallRNA_fam[which.max(Nunique)],
		                                                            forward_dir = forward_dir[which.max(Nunique)],
		                                                            smallrna_nomm = smallrna_nomm[which.max(Nunique)],
		                                                            mrna_length = mrna_length[which.max(Nunique)],
		                                                            cov = sum(cov),
		                                                            Nreads = sum(Nreads),
		                                                            mrna_aln = mrna_aln[which.max(Nunique)],
		                                                            mrna_aln2 = mrna_aln2[which.max(Nunique)],
		                                                            mrna_aln5 = mrna_aln5[which.max(Nunique)],
		                                                            mrna_aln10 = mrna_aln10[which.max(Nunique)],
		                                                            mrna_aln20 = mrna_aln20[which.max(Nunique)],
		                                                            seq.m = seq.m[which.max(Nunique)],
		                                                            Ndups=sum(Ndups),
		                                                            Nunique=sum(Nunique),
		                                                            noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                            noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)])
		                                                         ,by = .(chr.g,new_pos,seq.g.len,noncoding_clust) ]
		} else {
		  pair_occurence_sequence_filt_collap_clust <- input_tab[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                            flag.g = flag.g[which.max(Ndups)],
		                                                            strand = strand[which.max(Ndups)],
		                                                            smallRNA_fam = smallRNA_fam[which.max(Ndups)],
		                                                            forward_dir = forward_dir[which.max(Ndups)],
		                                                            smallrna_nomm = smallrna_nomm[which.max(Ndups)],
		                                                            mrna_length = mrna_length[which.max(Ndups)],
		                                                            cov = sum(cov),
		                                                            Nreads = sum(Nreads),
		                                                            mrna_aln = mrna_aln[which.max(Ndups)],
		                                                            mrna_aln2 = mrna_aln2[which.max(Ndups)],
		                                                            mrna_aln5 = mrna_aln5[which.max(Ndups)],
		                                                            mrna_aln10 = mrna_aln10[which.max(Ndups)],
		                                                            mrna_aln20 = mrna_aln20[which.max(Ndups)],
		                                                            seq.m = seq.m[which.max(Ndups)],
		                                                            Ndups=sum(Ndups),
		                                                            noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                            noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)])
		                                                         ,by = .(chr.g,new_pos,seq.g.len,noncoding_clust) ]
		}  
		pair_occurence_sequence_filt_collap_clust <- pair_occurence_sequence_filt_collap_clust[order(chr.g,new_pos)]
		names(pair_occurence_sequence_filt_collap_clust)[names(pair_occurence_sequence_filt_collap_clust) == 'new_pos'] <- 'pos.g'

		#collapse by target position and identified cluster

		#find all targets in a region of +- 20 and extend their sequence based on the longest union of positions  
		pair_occurence_sequence_filt_collap_clust[,new_start:=c(-21,pos.g[-length(pos.g)]), by=.(chr.g,noncoding_clust)]
		pair_occurence_sequence_filt_collap_clust[,new_start:=fifelse(pos.g <= new_start+20, FALSE, TRUE)]
		pair_occurence_sequence_filt_collap_clust[,new_pos:=pos.g]
		for(i in which(! pair_occurence_sequence_filt_collap_clust$new_start)){
		  if(i > 1)
		    set(pair_occurence_sequence_filt_collap_clust, i, "new_pos", pair_occurence_sequence_filt_collap_clust$new_pos[i-1L])
		}
		if (is_umi == TRUE){
		  pair_occurence_sequence_filt_collap_clust_2 <- pair_occurence_sequence_filt_collap_clust[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                                                              end.g = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),new_pos + seq.g.len],
		                                                                                              smallRNA_fam = smallRNA_fam[which.max(Nunique)],
		                                                                                              forward_dir = as.logical(max(forward_dir)),
		                                                                                              Ndups = sum(Ndups),
		                                                                                              Nunique = sum(Nunique), 
		                                                                                              flag.g = flag.g[which.max(Nunique)],
		                                                                                              strand = strand[which.max(Nunique)],
		                                                                                              smallrna_nomm = smallrna_nomm[which.max(Nunique)], 
		                                                                                              mrna_length = mrna_length[which.max(Nunique)],
		                                                                                              cov = sum(cov),
		                                                                                              Nreads = sum(Nreads),
		                                                                                              mrna_aln = mrna_aln[which.max(Nunique)], 
		                                                                                              mrna_aln2 = mrna_aln2[which.max(Nunique)], 
		                                                                                              mrna_aln5 = mrna_aln5[which.max(Nunique)], 
		                                                                                              mrna_aln10 = mrna_aln10[which.max(Nunique)],
		                                                                                              mrna_aln20 = mrna_aln20[which.max(Nunique)],
		                                                                                              real_smallRNA_seq = seq.m[which.max(Nunique)],
		                                                                                              noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                                              noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)]
		  ),
		  by=.(chr.g,new_pos,noncoding_clust)]
		} else {
		  pair_occurence_sequence_filt_collap_clust_2 <- pair_occurence_sequence_filt_collap_clust[,.(smallrna = paste0(unique(smallrna),collapse = "|"),
		                                                                                              end.g = .SD[order(-new_pos)][which.max(new_pos+seq.g.len),new_pos + seq.g.len],
		                                                                                              smallRNA_fam = smallRNA_fam[which.max(Ndups)],
		                                                                                              forward_dir = as.logical(max(forward_dir)),
		                                                                                              Ndups = sum(Ndups),
		                                                                                              flag.g = flag.g[which.max(Ndups)],
		                                                                                              strand = strand[which.max(Ndups)],
		                                                                                              smallrna_nomm = smallrna_nomm[which.max(Ndups)], 
		                                                                                              mrna_length = mrna_length[which.max(Ndups)],
		                                                                                              cov = sum(cov),
		                                                                                              Nreads = sum(Nreads),
		                                                                                              mrna_aln = mrna_aln[which.max(Ndups)], 
		                                                                                              mrna_aln2 = mrna_aln2[which.max(Ndups)], 
		                                                                                              mrna_aln5 = mrna_aln5[which.max(Ndups)], 
		                                                                                              mrna_aln10 = mrna_aln10[which.max(Ndups)],
		                                                                                              mrna_aln20 = mrna_aln20[which.max(Ndups)],
		                                                                                              real_smallRNA_seq = seq.m[which.max(Ndups)],
		                                                                                              noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
		                                                                                              noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)]
		  ),
		  by=.(chr.g,new_pos,noncoding_clust)]
	}
	pair_occurence_sequence_filt_collap_clust_2 <- pair_occurence_sequence_filt_collap_clust_2[order(chr.g,new_pos)]
	names(pair_occurence_sequence_filt_collap_clust_2)[names(pair_occurence_sequence_filt_collap_clust_2) == 'new_pos'] <- 'start.g'
	names(pair_occurence_sequence_filt_collap_clust_2)[names(pair_occurence_sequence_filt_collap_clust_2) == 'strand'] <- 'strand.g'
	return(pair_occurence_sequence_filt_collap_clust_2)
	}

	ambiguous <- pair_occurence_sequence_filt_collap_name[grepl("|",smallrna,fixed = T)]
	unique <- pair_occurence_sequence_filt_collap_name[!grepl("|",smallrna,fixed = T)]
	ambiguous_clust <- clustering(ambiguous)
	unique_clust <- clustering(unique)

	#####################################################
	merged_tab <- rbind(ambiguous_clust,unique_clust)
	pair_occurence_sequence_filt_collap <- merged_tab[order(chr.g,start.g)]
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'start.g'] <- 'new_pos'
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'end.g'] <- 'seq.g.len'
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'strand.g'] <- 'strand'
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'real_smallRNA_seq'] <- 'seq.m'
	pair_occurence_sequence_filt_collap[, seq.g.len := seq.g.len - new_pos]
	pair_occurence_sequence_filt_collap$noncoding_clust <- NULL
	merged_clust <- clustering(pair_occurence_sequence_filt_collap)
	pair_occurence_sequence_filt_collap <- merged_clust[order(chr.g,start.g)]


	#####################################################
	#to compute duplication ratio
	if (is_umi == TRUE){
		pair_occurence_sequence_filt_collap[,dups_perc := round(100-(Nunique/Ndups*100),1)]
	}
   
	#add genomic sequence of target alignment
	chr <- c(1:22,"MT","X","Y")
	pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[chr.g %in% chr]

	chromosome <- paste0("chr",pair_occurence_sequence_filt_collap$chr.g)
	chromosome <- gsub("chrMT","chrM",chromosome)

	ranges <- GRanges( chromosome,
	                 IRanges(start=pair_occurence_sequence_filt_collap$start.g, end = pair_occurence_sequence_filt_collap$end.g),
	                 strand=pair_occurence_sequence_filt_collap$strand.g )
	target_sequences <- as.data.table(Views(Hsapiens,ranges))
	#to retrieve original sequence from sequunced read
	pair_occurence_sequence_filt_collap[,seq.g := target_sequences$dna]



	##############################################################################################################################################################################################################################
	#####################################
	#check repeat_masker
	#add end position for repeatmasker overlap
	repeatmasker <- fread(repeatmasker_bed)
	names(repeatmasker) <- c("seqnames","start","end","repeatmasker","fam","type")
	repeatmasker <- repeatmasker[,.(seqnames, start, end, repeatmasker, fam)]
	setkeyv(repeatmasker, c("seqnames","start","end"))
	setkeyv(pair_occurence_sequence_filt_collap, c("chr.g","start.g","end.g"))

	overlapped_repeatmasker <- foverlaps(pair_occurence_sequence_filt_collap, repeatmasker, by.x=c("chr.g","start.g","end.g"), by.y=c("seqnames","start","end"))

	#check the overlap length of annotated repeatmasker
	overlapped_repeatmasker[,overlap:=apply(overlapped_repeatmasker[,c(2,3,6,9)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]

	#when we have more than 1 repeatmasker annotation at the same position, then take only those with longer overlap 
	res <- overlapped_repeatmasker[, .(paste0(unique(.SD[overlap == max(overlap),repeatmasker]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),fam]),collapse = ",")), by=.(chr.g,start.g,end.g)]
	res2 <- merge(res, overlapped_repeatmasker,by = c("chr.g","start.g","end.g"), all.y = T)
	if (is_umi == TRUE) {
		pair_occurence_sequence_filt_collap <- res2[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, Nunique, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type, noncodingRNA_seq, dups_perc, seq.g, V1, V2)]
	} else {
		pair_occurence_sequence_filt_collap <- res2[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type, noncodingRNA_seq, seq.g, V1, V2)]
	}
  	pair_occurence_sequence_filt_collap <- unique(pair_occurence_sequence_filt_collap)
  	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'V1'] <- 'repeatmasker'
  	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'V2'] <- 'repeat_fam'
	  
	#add gene annotation
	#########################
	annotate_by <- "gene_name"
	############load gtf file
	ref <- as.data.table(rtracklayer::import(transcripts))[type == "transcript", c("seqnames","start","end",annotate_by,"strand","gene_biotype"), with=F]

	setkeyv(ref, c("seqnames","start","end"))
	setkeyv(pair_occurence_sequence_filt_collap, c("chr.g","start.g","end.g"))

	# do overlap of genomic regions from alignment with gtf annotation file
	overlapped <- foverlaps(pair_occurence_sequence_filt_collap,ref,by.x=c("chr.g","start.g","end.g"),by.y=c("seqnames","start","end"),nomatch = 0)

	#check the overlap length of annotated genes
	overlapped[,overlap:=apply(overlapped[,c(2,3,7,8)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]

	res <- overlapped[, .(paste0(unique(.SD[overlap == max(overlap),gene_name]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),gene_biotype]),collapse = ",")), by=.(chr.g,start.g,end.g)]
	#add read name
	pair_occurence_sequence_filt_collap <- merge(res, pair_occurence_sequence_filt_collap, by = c("chr.g","start.g","end.g"), all.y = T)
	pair_occurence_sequence_filt_collap[,gene_name := V1]
	pair_occurence_sequence_filt_collap[,gene_biotype := V2]
	pair_occurence_sequence_filt_collap$V1 <- NULL
	pair_occurence_sequence_filt_collap$V2 <- NULL


	#add annotation of genomic feature location
	#########################
	annot_feat <- rtracklayer::import(features)

	annot_feat <- as.data.table(annot_feat)
	annot_feat <- annot_feat[,.(seqnames, start, end, type)]
	names(annot_feat) <- c("seqnames","start","end","loc")
	setkeyv(annot_feat, c("seqnames","start","end"))
	#add end position of target alignment
	setkeyv(pair_occurence_sequence_filt_collap, c("chr.g","start.g","end.g"))


	overlapped_feat2 <- foverlaps(pair_occurence_sequence_filt_collap, annot_feat, by.x=c("chr.g","start.g","end.g"), by.y=c("seqnames","start","end"))

	res_annot <- overlapped_feat2[, .(paste0(unique(loc),collapse = ",")), by=.(chr.g,start.g,end.g)]
	res_annot2 <- merge(res_annot, overlapped_feat2,by = c("chr.g","start.g","end.g"), all.y = T)
	if (is_umi == TRUE) {
		pair_occurence_sequence_filt_collap <- res_annot2[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, Nunique, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type,noncodingRNA_seq, dups_perc, seq.g, repeatmasker, repeat_fam, gene_name,gene_biotype, V1)]
  	} else {
    	pair_occurence_sequence_filt_collap <- res_annot2[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type,noncodingRNA_seq, seq.g, repeatmasker, repeat_fam, gene_name,gene_biotype, V1)]
	}
	pair_occurence_sequence_filt_collap <- unique(pair_occurence_sequence_filt_collap)

	#rename all added columns
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'V1'] <- 'feature'


	#########################################
	#fix smallrna names by real sequence

	fixed_smallRNA_names <- pair_occurence_sequence_filt_collap[,paste0(unique(unlist(strsplit(smallrna,"\\|"))),collapse = "|"),by = .(real_smallRNA_seq)]
	pair_occurence_sequence_filt_collap <- merge(pair_occurence_sequence_filt_collap,fixed_smallRNA_names, by ="real_smallRNA_seq")
	pair_occurence_sequence_filt_collap$smallrna <- NULL
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'V1'] <- 'smallrna'
	if (is_umi == TRUE) {
    	pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, Nunique, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type, noncodingRNA_seq, dups_perc, seq.g, repeatmasker, repeat_fam, gene_name, gene_biotype, feature)]
	} else {
    	pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[,.(chr.g, start.g, end.g, smallrna, smallRNA_fam, forward_dir, Ndups, flag.g, strand.g, smallrna_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, real_smallRNA_seq, noncodingRNA_type, noncodingRNA_seq, seq.g, repeatmasker, repeat_fam, gene_name, gene_biotype, feature)]
	}
	
	#######################################
	#update smallrna_nomm to small_nomm by checking all possible annotations
	pair_occurence_sequence_filt_collap[, smallrna_nomm := sapply(seq_along(real_smallRNA_seq),function(x) any(grepl(real_smallRNA_seq[x],wholedb[seq_name %in% unlist(tstrsplit(smallrna[x],"\\|"))]$sequence)))]
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'smallrna_nomm'] <- 'smallRNA_nomm'

  

	seqs <- sapply(seq_along(pair_occurence_sequence_filt_collap$smallrna), function(x) paste0(pair_occurence_sequence_filt_collap[x,seq.g],"&",if(pair_occurence_sequence_filt_collap[x,noncodingRNA_type == "miRNA" & !is.na(noncodingRNA_seq)]){pair_occurence_sequence_filt_collap[x,noncodingRNA_seq]}else{pair_occurence_sequence_filt_collap[x,real_smallRNA_seq]}))
	write.table(pair_occurence_sequence_filt_collap,paste0(c("hyb_pairs/",sample,"._before_cofold.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
	fa <- character(2 * length(seqs))
	fa[c(TRUE, FALSE)] = sprintf(">%s", seqs)
	fa[c(FALSE, TRUE)] = seqs
	writeLines(fa, paste0(c("hyb_pairs/cofold_pairs_collapsed/",sample,".filtered_collapsed_hybrid_sequence.fa"),collapse = ""))
	system(paste0(c("( cd hyb_pairs/cofold_pairs_collapsed/",sample,"/", " && RNAcofold --noLP --output-format='D' < ../",sample,".filtered_collapsed_hybrid_sequence.fa > ../",sample,".collapsed_cofold.csv )"),collapse = ""))
	
	pair_occurence_sequence_filt_collap[, A_perc := round(stri_count(real_smallRNA_seq, fixed = "A")/nchar(real_smallRNA_seq)*100,1)]
	pair_occurence_sequence_filt_collap[, C_perc := round(stri_count(real_smallRNA_seq, fixed = "C")/nchar(real_smallRNA_seq)*100,1)]
	pair_occurence_sequence_filt_collap[, G_perc := round(stri_count(real_smallRNA_seq, fixed = "G")/nchar(real_smallRNA_seq)*100,1)]
	pair_occurence_sequence_filt_collap[, T_perc := round(stri_count(real_smallRNA_seq, fixed = "T")/nchar(real_smallRNA_seq)*100,1)]


	cofold <- fread(paste0(c("hyb_pairs/cofold_pairs_collapsed/",sample,".collapsed_cofold.csv"),collapse = ""))
	cofold[,rna := tstrsplit(seq_id,"&")[1]]
	cofold[,smallrna := tstrsplit(seq_id,"&")[2]]
	pair_occurence_sequence_filt_collap<- cbind(pair_occurence_sequence_filt_collap,cofold[,.(mfe_struct,mfe)])
	pair_occurence_sequence_filt_collap[,ene :=5*(mfe-(-11))/(-5) ]
	pair_occurence_sequence_filt_collap[,ene := fifelse(ene >= 5,5,ene)]
	pair_occurence_sequence_filt_collap[,ene := fifelse(ene <= 0,0,ene)]
	pair_occurence_sequence_filt_collap[,noncoding_mirna_annot := grepl("hsa",smallrna)]
	pair_occurence_sequence_filt_collap[,noncoding_trna_annot := grepl("_tRNA",smallrna)]
	pair_occurence_sequence_filt_collap[,noncoding_rrna_annot := grepl("_rRNA",smallrna)]
	pair_occurence_sequence_filt_collap[,noncoding_snorna_annot := grepl("_snoRNA",smallrna)]
	pair_occurence_sequence_filt_collap[,noncoding_yrna_annot := grepl("_Y_RNA",smallrna)]
	pair_occurence_sequence_filt_collap[,noncoding_vaultrna_annot := grepl("_vaultRNA",smallrna)]

	pair_occurence_sequence_filt_collap[,noncoding_unique := sapply(seq_along(noncoding_mirna_annot), function(x) fifelse(sum(noncoding_mirna_annot[x],noncoding_trna_annot[x],noncoding_rrna_annot[x],noncoding_snorna_annot[x],noncoding_yrna_annot[x],noncoding_vaultrna_annot[x]) == 1,TRUE,FALSE))]

	#to fix noncoding types based on grepl
	pair_occurence_sequence_filt_collap[, noncodingRNA_type := gsub("^\\|","",paste0(fifelse(noncoding_mirna_annot,"|miRNA",""),fifelse(noncoding_rrna_annot,"|rRNA",""),fifelse(noncoding_snorna_annot,"|snoRNA",""),fifelse(noncoding_trna_annot,"|tRNA",""),fifelse(noncoding_yrna_annot,"|YRNA",""),fifelse(noncoding_vaultrna_annot,"|vaultRNA","")))]
	pair_occurence_sequence_filt_collap[, noncoding_simplified := tstrsplit(smallrna,"\\|")[1]]

	pair_occurence_sequence_filt_collap[, smallrna := sapply(seq_along(smallrna), function(x) paste0(unique(tstrsplit(smallrna[x],"\\|")),collapse = "|"))] 

	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'real_smallRNA_seq'] <- 'noncodingRNA_real_seq'
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'smallrna'] <- 'noncodingRNA'
	if (is_umi == TRUE) {
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, smallRNA_fam, forward_dir, Ndups, Nunique,
		                                                 dups_perc, repeatmasker, repeat_fam, smallRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
		                                                 mfe, ene, noncodingRNA_seq,noncoding_unique)] 
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[order(-Nunique)]
	} else {
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, smallRNA_fam, forward_dir, Ndups, 
		                                                 repeatmasker, repeat_fam, smallRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
		                                                 mfe, ene, noncodingRNA_seq,noncoding_unique)] 
		pair_occurence_sequence_filt_collap <- pair_occurence_sequence_filt_collap[order(-Ndups)]
	}
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'smallRNA_fam'] <- 'noncodingRNA_fam'
	names(pair_occurence_sequence_filt_collap)[names(pair_occurence_sequence_filt_collap) == 'smallRNA_nomm'] <- 'noncodingRNA_nomm'
	ambiguous_tab <- pair_occurence_sequence_filt_collap[noncoding_unique != T]
	unique_tab <- pair_occurence_sequence_filt_collap[noncoding_unique == T]
	ambiguous_tab$noncoding_unique <- NULL
	unique_tab$noncoding_unique <- NULL
	#this file contains all detected chimeric pairs deduplicated filtered and collapsed
	write.table(unique_tab,paste0(c("hyb_pairs/",sample,".hybrids_deduplicated_filtered_collapsed_unique.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
	write.table(ambiguous_tab,paste0(c("hyb_pairs/",sample,".hybrids_deduplicated_filtered_collapsed_ambiguous.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)

}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)