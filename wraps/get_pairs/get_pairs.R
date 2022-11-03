library(data.table)
library(seqinr)
library(miRBaseConverter)
library(Biostrings)
library(stringi)

run_all <- function(args){
	smallrna_bam <- args[1]
	sample <- args[2]
	rna_bam <- args[3]
	mod_rna_bed <- args[4]
    repeatmasker_bed <-args[5]
    smallrna_whole_bam <- args[6]
    transcripts <- args[7]
    features <- args[8]
    mirna_db <- args[9]
    is_umi <- args[10]
	smallrnabam <- fread(smallrna_bam)
	rnabam <- fread(rna_bam)
	rnabed <- fread(mod_rna_bed)
	repeatmasker <- fread(repeatmasker_bed)

	# necessary to determine position of miRNA and target
	smallrna_whole_bam <- fread(smallrna_whole_bam)

	#get read names from miRNA part, which aligned to the reverse strand and remve them from rnabam and smallrnabam
	reads_to_remove <- smallrnabam[V2 == 16 | V2 == 272]$V1
	rnabam <- rnabam[!V1 %in% reads_to_remove]
	smallrnabam <- smallrnabam[!V1 %in% reads_to_remove]

	# prepare rna genomic parts to annotate with gtf ref annotation
	#only kept the really aligned part of sequence for the annotation

	rnabam[,LSoftLen := gsub("^(\\d+)S.*","\\1",V6)]
	rnabam[,RSoftLen := gsub("^.*?(\\d+)S$","\\1",V6)]
	rnabam[,LSoftLen := as.numeric(LSoftLen)]
	rnabam[,RSoftLen := as.numeric(RSoftLen)]

	rnabam[is.na(rnabam)] <- 0
	rnabam[,seq.g_len := nchar(stri_sub(rnabam$V10, rnabam$LSoftLen+1 , (nchar(rnabam$V10)-rnabam$RSoftLen)))]
	comp_tab <- rnabam[,.(V3, V4, seq.g_len)]
	names(comp_tab) <- c("chr","start","aln")
	comp_tab[,end := start+aln]
	comp_tab[,name := rnabam[,V1]]
	comp_tab$aln <- NULL
	setkeyv(comp_tab, c("chr","start","end"))
	comp_tab[,chr:=as.character(chr)]

	#load annotation
	#########################
	############load gtf file
	ref <- as.data.table(rtracklayer::import(transcripts))[type == "transcript", c("seqnames","start","end","gene_name","strand"), with=F]
	setkeyv(ref, c("seqnames","start","end"))

	# do overlap of genomic regions from alignment with gtf annotation file
	#overlapped <- comp_tab[ref,on=c("chr==seqnames"),allow.cartesian=T][start+1<=i.end & end>=i.start]
	overlapped <- foverlaps(comp_tab,ref,by.x=c("chr","start","end"),by.y=c("seqnames","start","end"),nomatch = 0)

	overlapped[,c("start","end","i.start","i.end") := list(i.start, i.end, start, end)]

	#apply(overlapped[,c(2,3,6,7)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)

	#check the overlap length of annotated genes
	overlapped[,overlap:=apply(overlapped[,.(start, end, i.start, i.end)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]

	res <- overlapped[, .(paste0(unique(.SD[overlap == max(overlap),gene_name]),collapse = ",")), by=.(chr,start,end)]
	#resolve targets with overlapping annotation
	res[,gene_name:=gsub("(A[A-Z][0-9]{6}(\\.[0-9])?,|,A[A-Z][0-9]{6}(\\.[0-9])?$)","",V1)]
	res$V1 <- NULL
	#add read name
	res <- merge(res, comp_tab, by = c("chr","start","end"), all.y = T)

	#########################
	# merge together rna and smallrna bam with anotation of rna parts
	annot_tab_new <- res[,-3]
	#add gene_name annotation to rnapart bam
	#keep read name, flag, chromosome, position, mapping quality,cigar and sequence
	rnabam_gene_new <- merge(rnabam[,.(V1, V2, V3, V4, V5, V6, V10)],annot_tab_new[,c(3,4)], by.x = "V1",by.y = "name")

	#add strand sign to each chimera
	#the reads with flag 0 are mapped to the forward strand and the reads with flag 16 are mapped to the reverse strand.
	check_strandness <- function(x){
	  if (intToBits(x)[5] == 1){
	    return("-")
	  }else{
	    return("+")
	  }
	}

	rnabam_gene_new[,strand := sapply(V2,function(x) check_strandness(x))]

	# add to RNA part bam also the information about miRNA mapping - read name, smallrna, position, mapping quality,cigar and sequence
	candidates_gene_new <- merge(rnabam_gene_new,smallrnabam[,.(V1, V3, V4, V5, V6, V10)],by = "V1")  
	colnames(candidates_gene_new) <- c("read_name","flag.g","chr.g","pos.g","qual.g","cigar.g","seq.g","gene_name","strand","smallrna","pos.m","qual.m","cigar.m","seq.m")

	#create table of pairs based just on gene and smallrna names and their counts
	pair_occurence_gene <- candidates_gene_new[,.N,by=.(gene_name,smallrna)]

	# add the sequence of whole read 
	#check direction of original read
	#forward_dir is true when the direction was mRNA-miRNA
	smallrna_whole_bam[,LSoftLen := as.numeric(gsub("^(\\d+)S.*","\\1",V6))]
	smallrna_whole_bam[,RSoftLen := as.numeric(gsub("^.*?(\\d+)S$","\\1",V6))]
	smallrna_whole_bam[is.na(smallrna_whole_bam)] <- 0
	#take only reads, where one of softclip site is longer than 15 bp and matched part for miRNA need to be longer than 15bp
	smallrna_whole_bam[,forward_dir := fifelse(LSoftLen >= RSoftLen,TRUE,FALSE)]
	candidates_gene_new <- merge(candidates_gene_new,smallrna_whole_bam[,.(V1, forward_dir)],by.x = "read_name",by.y = "V1")


	if (is_umi == TRUE) {
		# extract UMI for deduplication
		candidates_gene_new$UMI <- unlist(strsplit(candidates_gene_new$read_name, "_"))[seq(2,dim(candidates_gene_new)[1]*2,2)]
		candidates_gene_new <- candidates_gene_new[order(candidates_gene_new$UMI)]

		#number of different UMIs at same aln position
		collap <- candidates_gene_new[,length(unique(UMI)),by=.(chr.g,flag.g,pos.g,cigar.g,gene_name,seq.g,seq.m,smallrna,strand,forward_dir)]
		#here we get the total number of reads aligned to the same aln position, regardless the UMI
		#PLEASE CHECK HERE if with new smallrnawhole_bam, number of collapsed reads here remain the same with and without forward_dir column
		pair_occurence_sequence <- candidates_gene_new[,.N,by=.(chr.g,flag.g,pos.g,cigar.g,gene_name,seq.g,seq.m,smallrna,strand,forward_dir)]


		# add the number of different UMIs at the same aln position
		  pair_occurence_sequence <- merge(pair_occurence_sequence,collap,by=c("chr.g", "flag.g", "pos.g", "cigar.g", "gene_name", "seq.g", "seq.m", "smallrna", "strand", "forward_dir"))

		names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'N'] <- 'Ndups'
		names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V1'] <- 'Nunique'

		#add percentage of duplicates
		pair_occurence_sequence[,dups_perc:=round(100 - (Nunique/Ndups*100),1)]
	  
	} else {
		pair_occurence_sequence <- candidates_gene_new[,.N,by=.(chr.g,flag.g,pos.g,cigar.g,gene_name,seq.g,seq.m,smallrna,strand,forward_dir)]
		# add the number of different UMIs at the same aln position
		#pair_occurence_sequence <- merge(pair_occurence_sequence,b,by=colnames(pair_occurence_sequence[,1:10]))

		names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'N'] <- 'Ndups'
	}


	####################################
	#####add mirna family
	mirna_accessions <- miRNA_NameToAccession(pair_occurence_sequence$smallrna,version = "v22")
	mirna_fam <- checkMiRNAFamily(mirna_accessions$Accession)

	#add mirna family annotation column
	pair_occurence_sequence <- cbind(pair_occurence_sequence,mirna_fam$Family) 


	####################################
	# take only mapped part of target sequence and remove softclipped parts
	pair_occurence_sequence[,LSoftLen := as.numeric(gsub("^(\\d+)S.*","\\1",cigar.g))]
	pair_occurence_sequence[,RSoftLen := as.numeric(gsub("^.*?(\\d+)S$","\\1",cigar.g))]
	pair_occurence_sequence[is.na(pair_occurence_sequence)] <- 0
	#do reverse complement on mRNA seqeunces, where there is flag 16
	#change seq.g to reverse complement in a case, where we see flag 16 
	#this construct is terribly slow, as the function is called many times
	#pair_occurence_sequence[flag.g == 16 , seq.g := sapply(seq.g,function(x) as.character(reverseComplement(DNAStringSet(c(x)))[[1]]))]
	pair_occurence_sequence[flag.g == 16 , seq.g := as.character(reverseComplement(DNAStringSet(seq.g)))]
	#switch softclip parts, to have them properly cutted
	pair_occurence_sequence[flag.g == 16,c("LSoftLen", "RSoftLen") := .(RSoftLen, LSoftLen)]
	pair_occurence_sequence[,seq.g := stri_sub(pair_occurence_sequence$seq.g, pair_occurence_sequence$LSoftLen+1 , (nchar(pair_occurence_sequence$seq.g)-pair_occurence_sequence$RSoftLen))]
	pair_occurence_sequence$LSoftLen <- NULL
	pair_occurence_sequence$RSoftLen <- NULL

	####################################
	#join sequences in a direction either mRNA-miRNA or miRNA-mRNA
	pair_occurence_sequence[,read_seq := fifelse(forward_dir == TRUE,paste0(seq.g,seq.m),paste0(seq.m,seq.g))]

	#####################################
	#check repeat_masker
	#add end position for repeatmasker overlap
	pair_occurence_sequence[,end := pos.g + nchar(seq.g)]
	names(repeatmasker) <- c("seqnames","start","end","repeatmasker","fam","type")
	repeatmasker <- repeatmasker[,.(seqnames, start, end, repeatmasker, fam)]
	setkeyv(repeatmasker, c("seqnames","start","end"))
	setkeyv(pair_occurence_sequence, c("chr.g","pos.g","end"))
	overlapped_repeatmasker <- foverlaps(pair_occurence_sequence, repeatmasker, by.x=c("chr.g","pos.g","end"), by.y=c("seqnames","start","end"))

	
	if (is_umi == TRUE){
		#check the overlap length of annotated repeatmasker
 		overlapped_repeatmasker[,overlap:=apply(overlapped_repeatmasker[,.(start, end, pos.g, i.end)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]

		#when we have more than 1 repeatmasker annotation at the same position, then take only those with longer overlap 		
		res <- overlapped_repeatmasker[, .(paste0(unique(.SD[overlap == max(overlap),repeatmasker]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),fam]),collapse = ",")), by=.(chr.g,pos.g,i.end)]
		res2 <- merge(res, overlapped_repeatmasker,by = c("chr.g","pos.g","i.end"), all.y = T)
		pair_occurence_sequence <- res2[,c(1,2,10:22,4,5)]
	} else{
		#check the overlap length of annotated repeatmasker
		overlapped_repeatmasker[,overlap:=apply(overlapped_repeatmasker[,.(start, end, pos.g, i.end)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]

		#when we have more than 1 repeatmasker annotation at the same position, then take only those with longer overlap 
		res <- overlapped_repeatmasker[, .(paste0(unique(.SD[overlap == max(overlap),repeatmasker]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),fam]),collapse = ",")), by=.(chr.g,pos.g,i.end)]
		res2 <- merge(res, overlapped_repeatmasker,by = c("chr.g","pos.g","i.end"), all.y = T)
		pair_occurence_sequence <- res2[,.(chr.g, pos.g, flag.g, cigar.g, gene_name, seq.g, seq.m, smallrna, strand, forward_dir, Ndups, V2.y, read_seq, V1, V2.x)]
	}

	
	pair_occurence_sequence <- unique(pair_occurence_sequence)

	#add annotation of genomic feature location
	annot_feat <- rtracklayer::import(features)

	annot_feat <- as.data.table(annot_feat)
	annot_feat <- annot_feat[,.(seqnames, start, end, type)]
	names(annot_feat) <- c("seqnames","start","end","loc")
	setkeyv(annot_feat, c("seqnames","start","end"))
	#add end position of target alignment
	pair_occurence_sequence[,end := pos.g + nchar(seq.g)] 
	setkeyv(pair_occurence_sequence, c("chr.g","pos.g","end"))


	overlapped_feat2 <- foverlaps(pair_occurence_sequence, annot_feat, by.x=c("chr.g","pos.g","end"), by.y=c("seqnames","start","end"))

	res_annot <- overlapped_feat2[, .(paste0(unique(loc),collapse = ",")), by=.(chr.g,pos.g,i.end)]
	res_annot2 <- merge(res_annot, overlapped_feat2,by = c("chr.g","pos.g","i.end"), all.y = T)
	if (is_umi == TRUE) {
		pair_occurence_sequence <- res_annot2[,.(chr.g, pos.g, flag.g, cigar.g, gene_name, seq.g, seq.m, smallrna, strand, forward_dir, Ndups, Nunique, dups_perc, V2.y, read_seq, V1.y, V2.x, V1.x)]
	} else {
		pair_occurence_sequence <- res_annot2[,.(chr.g, pos.g, flag.g, cigar.g, gene_name, seq.g, seq.m, smallrna, strand, forward_dir, Ndups, V2.y, read_seq, V1.y, V2.x, V1.x)]
	}		
	pair_occurence_sequence <- unique(pair_occurence_sequence)

	#rename all added columns
	names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V2.y'] <- 'smallRNA_fam'
	names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V1.y'] <- 'repeatmasker'
	names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V2.x'] <- 'fam'
	names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V1.x'] <- 'feature'

	#######
	#######
	#add filtering for 0 mismatch miRNA alignment
	mirnadb <- read.fasta(mirna_db,as.string = TRUE)
	seq_name <- names(mirnadb)
	sequence <- paste(mirnadb)
	mirnaDB <- data.table(seq_name, sequence)
	mirnaDB[, sequence := toupper(sequence)]


	pair_occurence_sequence[, smallrna_nomm := sapply(seq.m,function(x) any(grepl(x,mirnaDB$sequence)))]

	pair_occurence_sequence[, mrna_length := nchar(seq.g) >= 20]

	#####################################
	#####################################
	#check genomic mapping
	names(rnabed) <- c("seqnames","start","end","cov")
	setkeyv(rnabed, c("seqnames","start","end"))
	rnabed <- rnabed[-which(rnabed$end < rnabed$start)]
	pair_occurence_sequence[,end := pos.g + nchar(seq.g)] 

	overlapped_rna_aln <- foverlaps(pair_occurence_sequence, rnabed, by.x=c("chr.g","pos.g","end"), by.y=c("seqnames","start","end"))
	pair_occurence_sequence$end <- NULL

	res_aln <- overlapped_rna_aln[, .(max(cov)), by=.(chr.g,pos.g,i.end)]
	res_aln2 <- merge(res_aln, overlapped_rna_aln,by = c("chr.g","pos.g","i.end"), all.y = T)

	if (is_umi == TRUE) {
		pair_occurence_sequence <- res_aln2[,.(chr.g, pos.g, flag.g, cigar.g, gene_name, seq.g, seq.m, smallrna, strand, forward_dir, Ndups, Nunique, dups_perc, smallRNA_fam,read_seq, repeatmasker, fam, feature, mrna_length, smallrna_nomm, V1)]
	} else {
		pair_occurence_sequence <- res_aln2[,.(chr.g, pos.g, flag.g, cigar.g, gene_name, seq.g, seq.m, smallrna, strand, forward_dir, Ndups, smallRNA_fam,read_seq, repeatmasker, fam, feature, mrna_length, smallrna_nomm, V1)]
	}	
	
	pair_occurence_sequence <- unique(pair_occurence_sequence)

	#rename all added columns
	names(pair_occurence_sequence)[names(pair_occurence_sequence) == 'V1'] <- 'cov'

	# define all chimeras, which are supported by at least one read on genomic location by simple reads
	pair_occurence_sequence[,mrna_aln := !is.na(cov)]

	#genomic mapping with at least 2 coverage coverage
	pair_occurence_sequence[, Nreads := round(cov/min(pair_occurence_sequence[!is.na(cov),]$cov),0),]


	pair_occurence_sequence[,mrna_aln2 := !is.na(cov) & Nreads >= 2]
	pair_occurence_sequence[,mrna_aln5 := !is.na(cov) & Nreads >= 5]
	pair_occurence_sequence[,mrna_aln10 := !is.na(cov) & Nreads >= 10]
	pair_occurence_sequence[,mrna_aln20 := !is.na(cov) & Nreads >= 20]

	if (is_umi == TRUE){
		pair_occurence_sequence_sort <- pair_occurence_sequence[order(-Nunique)]
	} else {
		pair_occurence_sequence_sort <- pair_occurence_sequence[order(-Ndups)]
	}
	
	seqs <- paste0(pair_occurence_sequence[,seq.g],"&",pair_occurence_sequence[,seq.m])
	write.table(pair_occurence_sequence,paste0(c("hyb_pairs/",sample,".sequence_hyb_pairs_deduplicated_before_cofold.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
	#this file contains all detected chimeric pairs without any deduplication or filtering and collapsing
	write.table(candidates_gene_new,paste0(c("hyb_pairs/",sample,".all_hyb_pairs.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
	fa <- character(2 * length(seqs))
	fa[c(TRUE, FALSE)] = sprintf(">%s", seqs)
	fa[c(FALSE, TRUE)] = seqs
	writeLines(fa, paste0(c("hyb_pairs/cofold_pairs/",sample,".all_hybrid_sequence.fa"),collapse = ""))
	system(paste0(c("( cd hyb_pairs/cofold_pairs/",sample,"/", " && RNAcofold --noLP --output-format='D' < ../",sample,".all_hybrid_sequence.fa > ../",sample,".all_hybrid_cofold.csv )"),collapse = ""))
	#add relative abundance
	pair_occurence_sequence[, A_perc := round(stri_count(seq.m, fixed = "A")/nchar(seq.m)*100,1)]
	pair_occurence_sequence[, C_perc := round(stri_count(seq.m, fixed = "C")/nchar(seq.m)*100,1)]
	pair_occurence_sequence[, G_perc := round(stri_count(seq.m, fixed = "G")/nchar(seq.m)*100,1)]
	pair_occurence_sequence[, T_perc := round(stri_count(seq.m, fixed = "T")/nchar(seq.m)*100,1)]

	cofold <- fread(paste0(c("hyb_pairs/cofold_pairs/",sample,".all_hybrid_cofold.csv"),collapse = ""))
	cofold[,rna := tstrsplit(seq_id,"&")[1]]
	cofold[,smallrna := tstrsplit(seq_id,"&")[2]]
	pair_occurence_sequence<- cbind(pair_occurence_sequence,cofold[,.(mfe_struct,mfe)])
	# convert to the energy numbers according their formula
	pair_occurence_sequence[,ene :=5*(mfe-(-11))/(-5) ]
	pair_occurence_sequence[,ene := fifelse(ene >= 5,5,ene)]
	pair_occurence_sequence[,ene := fifelse(ene <= 0,0,ene)]
	if (is_umi == TRUE){
		pair_occurence_sequence <- pair_occurence_sequence[order(-Nunique)]
	} else {
		pair_occurence_sequence <- pair_occurence_sequence[order(-Ndups)]
	}

	#this file contains all detected chimeric pairs deduplicated but without filtering and collapsing
	write.table(pair_occurence_sequence,paste0(c("hyb_pairs/",sample,".hyb_pairs_deduplicated.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)