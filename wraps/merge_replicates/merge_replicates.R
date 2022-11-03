library(data.table)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringi)


run_all <- function(args){
	args <- as.list(args) # making the list
	repeatmasker_bed <-args[1]
	transcripts <- args[2]
	features <- args[3]
	is_umi <- args[4]
	is_ambiguous <- args[5]


	for (i in 6:length(args)){
    sample <- fread(unlist(args[i]))
    sample[,rep := paste("rep",i-5,sep = "")]
    if (!exists("all_reps_tab")){
      all_reps_tab <- copy(sample)
    } else {
      all_reps_tab <- rbind(all_reps_tab,sample)
    }
  } 
  
  #all_reps <- copy(all_reps_mirna)
  repeatmasker <- fread(unlist(repeatmasker_bed))
  ref <- as.data.table(rtracklayer::import(unlist(transcripts)))[type == "transcript", c("seqnames","start","end","gene_name","strand","gene_biotype"), with=F]
  annot_feat <- rtracklayer::import(unlist(features))
  
  merge_and_collapse <- function(all_reps,repeatmasker,ref,annot_feat,type){
    ##############################
    #
    all_reps <- all_reps[order(chr.g,noncodingRNA,start.g)]
    all_reps[,new_start:=c(-21,start.g[-length(start.g)]), by=.(chr.g,noncodingRNA)]
    all_reps[,new_start:=fifelse(start.g <= new_start+20, FALSE, TRUE)]
    all_reps[,new_pos:=start.g]
    for(i in which(! all_reps$new_start)){
      if(i > 1)
        set(all_reps, i, "new_pos", all_reps$new_pos[i-1L])
    }
    
    #to check positions
    if (is_umi == TRUE) {
      all_reps_collap <- all_reps[,.(end.g = max(end.g), 
                                     Ndups=sum(Ndups),
                                     Nunique=sum(Nunique),
                                     gene_name=gene_name[which.max(Nunique)],
                                     repeatmasker=repeatmasker[which.max(Nunique)],
                                     repeat_fam=repeat_fam[which.max(Nunique)],
                                     feature=feature[which.max(Nunique)],
                                     noncodingRNA_nomm=noncodingRNA_nomm[which.max(Nunique)],
                                     mrna_length=mrna_length[which.max(Nunique)],
                                     cov=cov[which.max(Nunique)],
                                     Nreads=Nreads[which.max(Nunique)],
                                     mrna_aln=mrna_aln[which.max(Nunique)],
                                     mrna_aln2= mrna_aln2[which.max(Nunique)],
                                     mrna_aln5= mrna_aln5[which.max(Nunique)],
                                     mrna_aln10= mrna_aln10[which.max(Nunique)],
                                     mrna_aln20= mrna_aln20[which.max(Nunique)],
                                     rep = paste0(unique(rep),collapse = ","),
                                     noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Nunique)]
                                     
      ),
      by=.(chr.g,new_pos,flag.g,strand.g,noncodingRNA,noncodingRNA_seq,noncodingRNA_type,noncodingRNA_fam,forward_dir)]
    } else {
      all_reps_collap <- all_reps[,.(end.g = max(end.g),
                                     Ndups=sum(Ndups),
                                     gene_name=gene_name[which.max(Ndups)],
                                     repeatmasker=repeatmasker[which.max(Ndups)],
                                     repeat_fam=repeat_fam[which.max(Ndups)],
                                     feature=feature[which.max(Ndups)],
                                     noncodingRNA_nomm=noncodingRNA_nomm[which.max(Ndups)],
                                     mrna_length=mrna_length[which.max(Ndups)],
                                     cov=cov[which.max(Ndups)],
                                     Nreads=Nreads[which.max(Ndups)],
                                     mrna_aln=mrna_aln[which.max(Ndups)],
                                     mrna_aln2= mrna_aln2[which.max(Ndups)],
                                     mrna_aln5= mrna_aln5[which.max(Ndups)],
                                     mrna_aln10= mrna_aln10[which.max(Ndups)],
                                     mrna_aln20= mrna_aln20[which.max(Ndups)],
                                     rep = paste0(unique(rep),collapse = ","),
                                     noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Ndups)]
                                     
      ),
      by=.(chr.g,new_pos,flag.g,strand.g,noncodingRNA,noncodingRNA_seq,noncodingRNA_type,noncodingRNA_fam,forward_dir)]
    }
    
    
    
    #perform sequence similarity clustering
    distance_noncoding <- adist(all_reps_collap$noncodingRNA_real_seq,all_reps_collap$noncodingRNA_real_seq)
    mean_len_mat_noncoding <- as.matrix(proxy::dist(nchar(all_reps_collap$noncodingRNA_real_seq),method = function(x,y) (x+y)/2))
    diag(mean_len_mat_noncoding) <- nchar(all_reps_collap$noncodingRNA_real_seq)
    all_reps_collap[,noncoding_clust := cutree(hclust(as.dist(distance_noncoding/mean_len_mat_noncoding)),h = 0.3)]
    #collapse by identified clusters
    if (is_umi == TRUE) {
      pair_occurence_sequence_filt_collap_clust <- all_reps_collap[,.(noncodingRNA = paste0(unique(noncodingRNA),collapse = "|"),
                                                                      flag.g = flag.g[which.max(Nunique)],
                                                                      strand.g = strand.g[which.max(Nunique)],
                                                                      noncodingRNA_fam = noncodingRNA_fam[which.max(Nunique)],
                                                                      forward_dir = forward_dir[which.max(Nunique)],
                                                                      noncodingRNA_nomm = noncodingRNA_nomm[which.max(Nunique)],
                                                                      mrna_length = mrna_length[which.max(Nunique)],
                                                                      cov = cov[which.max(Nunique)],
                                                                      Nreads = Nreads[which.max(Nunique)],
                                                                      mrna_aln = mrna_aln[which.max(Nunique)],
                                                                      mrna_aln2 = mrna_aln2[which.max(Nunique)],
                                                                      mrna_aln5 = mrna_aln5[which.max(Nunique)],
                                                                      mrna_aln10 = mrna_aln10[which.max(Nunique)],
                                                                      mrna_aln20 = mrna_aln20[which.max(Nunique)],
                                                                      Ndups=sum(Ndups),
                                                                      Nunique=sum(Nunique),
                                                                      gene_name=gene_name[which.max(Nunique)],
                                                                      repeatmasker=repeatmasker[which.max(Nunique)],
                                                                      repeat_fam=repeat_fam[which.max(Nunique)],
                                                                      feature=feature[which.max(Nunique)],
                                                                      noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
                                                                      noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)],
                                                                      noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Nunique)],
                                                                      rep = paste0(unique(rep),collapse = ","))
                                                                   ,by = .(chr.g,new_pos,end.g,noncoding_clust) ]
    } else {
      pair_occurence_sequence_filt_collap_clust <- all_reps_collap[,.(noncodingRNA = paste0(unique(noncodingRNA),collapse = "|"),
                                                                      flag.g = flag.g[which.max(Ndups)],
                                                                      strand.g = strand.g[which.max(Ndups)],
                                                                      noncodingRNA_fam = noncodingRNA_fam[which.max(Ndups)],
                                                                      forward_dir = forward_dir[which.max(Ndups)],
                                                                      noncodingRNA_nomm = noncodingRNA_nomm[which.max(Ndups)],
                                                                      mrna_length = mrna_length[which.max(Ndups)],
                                                                      cov = cov[which.max(Ndups)],
                                                                      Nreads = Nreads[which.max(Ndups)],
                                                                      mrna_aln = mrna_aln[which.max(Ndups)],
                                                                      mrna_aln2 = mrna_aln2[which.max(Ndups)],
                                                                      mrna_aln5 = mrna_aln5[which.max(Ndups)],
                                                                      mrna_aln10 = mrna_aln10[which.max(Ndups)],
                                                                      mrna_aln20 = mrna_aln20[which.max(Ndups)],
                                                                      Ndups=sum(Ndups),
                                                                      gene_name=gene_name[which.max(Ndups)],
                                                                      repeatmasker=repeatmasker[which.max(Ndups)],
                                                                      repeat_fam=repeat_fam[which.max(Ndups)],
                                                                      feature=feature[which.max(Ndups)],
                                                                      noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
                                                                      noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)],
                                                                      noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Ndups)],
                                                                      rep = paste0(unique(rep),collapse = ","))
                                                                   ,by = .(chr.g,new_pos,end.g,noncoding_clust) ]
    }  
    
    names(pair_occurence_sequence_filt_collap_clust)[names(pair_occurence_sequence_filt_collap_clust) == 'new_pos'] <- 'pos.g'
    
    pair_occurence_sequence_filt_collap_clust <- pair_occurence_sequence_filt_collap_clust[order(chr.g,noncoding_clust,pos.g)]
    pair_occurence_sequence_filt_collap_clust[,new_start:=c(-21,pos.g[-length(pos.g)]), by=.(chr.g,noncoding_clust)]
    pair_occurence_sequence_filt_collap_clust[,new_start:=fifelse(pos.g <= new_start+20, FALSE, TRUE)]
    pair_occurence_sequence_filt_collap_clust[,new_pos:=pos.g]
    for(i in which(! pair_occurence_sequence_filt_collap_clust$new_start)){
      if(i > 1)
        set(pair_occurence_sequence_filt_collap_clust, i, "new_pos", pair_occurence_sequence_filt_collap_clust$new_pos[i-1L])
    }
    if (is_umi == TRUE) {
      pair_occurence_sequence_filt_collap_clust_2 <- pair_occurence_sequence_filt_collap_clust[,.(noncodingRNA = paste0(unique(noncodingRNA),collapse = "|"),
                                                                                                  end.g = max(end.g),
                                                                                                  noncodingRNA_fam = noncodingRNA_fam[which.max(Nunique)],
                                                                                                  forward_dir = as.logical(max(forward_dir)),
                                                                                                  Ndups = sum(Ndups),
                                                                                                  Nunique = sum(Nunique),
                                                                                                  flag.g = flag.g[which.max(Nunique)],
                                                                                                  strand.g = strand.g[which.max(Nunique)],
                                                                                                  noncodingRNA_nomm = noncodingRNA_nomm[which.max(Nunique)],
                                                                                                  mrna_length = mrna_length[which.max(Nunique)],
                                                                                                  cov = cov[which.max(Nunique)],
                                                                                                  Nreads = Nreads[which.max(Nunique)],
                                                                                                  mrna_aln = mrna_aln[which.max(Nunique)],
                                                                                                  mrna_aln2 = mrna_aln2[which.max(Nunique)],
                                                                                                  mrna_aln5 = mrna_aln5[which.max(Nunique)],
                                                                                                  mrna_aln10 = mrna_aln10[which.max(Nunique)],
                                                                                                  mrna_aln20 = mrna_aln20[which.max(Nunique)],
                                                                                                  noncodingRNA_seq = noncodingRNA_seq[which.max(Nunique)],
                                                                                                  noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
                                                                                                  gene_name=gene_name[which.max(Nunique)],
                                                                                                  repeatmasker=repeatmasker[which.max(Nunique)],
                                                                                                  repeat_fam=repeat_fam[which.max(Nunique)],
                                                                                                  feature=feature[which.max(Nunique)],
                                                                                                  noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Nunique)],
                                                                                                  rep = paste0(unique(rep),collapse = ",")
      ),
      by=.(chr.g,new_pos,noncoding_clust)]
      #to compute duplication ratio
      pair_occurence_sequence_filt_collap_clust_2[,dups_perc := round(100-(Nunique/Ndups*100),1)]
    } else {
      pair_occurence_sequence_filt_collap_clust_2 <- pair_occurence_sequence_filt_collap_clust[,.(noncodingRNA = paste0(unique(noncodingRNA),collapse = "|"),
                                                                                                  end.g = max(end.g),
                                                                                                  noncodingRNA_fam = noncodingRNA_fam[which.max(Ndups)],
                                                                                                  forward_dir = as.logical(max(forward_dir)),
                                                                                                  Ndups = sum(Ndups),
                                                                                                  flag.g = flag.g[which.max(Ndups)],
                                                                                                  strand.g = strand.g[which.max(Ndups)],
                                                                                                  noncodingRNA_nomm = noncodingRNA_nomm[which.max(Ndups)],
                                                                                                  mrna_length = mrna_length[which.max(Ndups)],
                                                                                                  cov = cov[which.max(Ndups)],
                                                                                                  Nreads = Nreads[which.max(Ndups)],
                                                                                                  mrna_aln = mrna_aln[which.max(Ndups)],
                                                                                                  mrna_aln2 = mrna_aln2[which.max(Ndups)],
                                                                                                  mrna_aln5 = mrna_aln5[which.max(Ndups)],
                                                                                                  mrna_aln10 = mrna_aln10[which.max(Ndups)],
                                                                                                  mrna_aln20 = mrna_aln20[which.max(Ndups)],
                                                                                                  noncodingRNA_seq = noncodingRNA_seq[which.max(Ndups)],
                                                                                                  noncodingRNA_type = paste0(unique(noncodingRNA_type),collapse = "|"),
                                                                                                  gene_name=gene_name[which.max(Ndups)],
                                                                                                  repeatmasker=repeatmasker[which.max(Ndups)],
                                                                                                  repeat_fam=repeat_fam[which.max(Ndups)],
                                                                                                  feature=feature[which.max(Ndups)],
                                                                                                  noncodingRNA_real_seq = noncodingRNA_real_seq[which.max(Ndups)],
                                                                                                  rep = paste0(unique(rep),collapse = ",")
      ),
      by=.(chr.g,new_pos,noncoding_clust)]
    }
    all_reps_collap_2 <- pair_occurence_sequence_filt_collap_clust_2[order(chr.g,new_pos)]
    
    names(all_reps_collap_2)[names(all_reps_collap_2) == 'new_pos'] <- 'start.g'
    
    chromosome <- paste0("chr",all_reps_collap_2$chr.g)
    chromosome <- gsub("chrMT","chrM",chromosome)
    
    ranges <- GRanges( chromosome,
                       IRanges(start=all_reps_collap_2$start.g, end = all_reps_collap_2$end.g),
                       strand=all_reps_collap_2$strand.g )
    target_sequences <- as.data.table(Views(Hsapiens,ranges))
    #here in target sequences we got reverse complement, but need only complement
    target_sequences[strand == "-", dna_rev := sapply(dna,function(x) as.character(reverse(DNAStringSet(c(x)))))]
    all_reps_collap_2[,seq.g := target_sequences$dna]
    all_reps <- all_reps[order(chr.g,noncodingRNA,start.g)]

    ##############################################################################################################################################################################################################################
    #####################################
    # fix_annotations 
    #check repeat_masker
    #add end position for repeatmasker overlap
    hybrids_tab <- copy(all_reps)
    names(repeatmasker) <- c("seqnames","start","end","repeatmasker","fam","type")
    repeatmasker <- repeatmasker[,.(seqnames,start,end,repeatmasker,fam)]
    setkeyv(repeatmasker, c("seqnames","start","end"))
    setkeyv(hybrids_tab, c("chr.g","start.g","end.g"))
    
    overlapped_repeatmasker <- foverlaps(hybrids_tab, repeatmasker, by.x=c("chr.g","start.g","end.g"), by.y=c("seqnames","start","end"))
    
    #check the overlap length of annotated repeatmasker
    overlapped_repeatmasker[,overlap:=apply(overlapped_repeatmasker[,.(start, end, start.g, end.g)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]
    
    #when we have more than 1 repeatmasker annotation at the same position, then take only those with longer overlap
    res <- overlapped_repeatmasker[, .(paste0(unique(.SD[overlap == max(overlap),repeatmasker]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),fam]),collapse = ",")), by=.(chr.g,start.g,end.g)]
    
    hybrids_tab <- merge(res, overlapped_repeatmasker,by = c("chr.g","start.g","end.g"), all.y = T)
    hybrids_tab[, i.repeatmasker := V1]
    hybrids_tab$V1 <- NULL
    hybrids_tab$repeatmasker <- NULL
    names(hybrids_tab)[names(hybrids_tab) == 'i.repeatmasker'] <- 'repeatmasker'
    hybrids_tab[, repeat_fam := V2]
    hybrids_tab$V2 <- NULL
    hybrids_tab$fam <- NULL
    hybrids_tab$start <- NULL
    hybrids_tab$end <- NULL
    hybrids_tab$overlap <- NULL
    hybrids_tab <- unique(hybrids_tab)
    #add gene annotation
    #########################
    ############load gtf file
    setkeyv(ref, c("seqnames","start","end"))
    setkeyv(hybrids_tab, c("chr.g","start.g","end.g"))
    
    # do overlap of genomic regions from alignment with gtf annotation file
    overlapped <- foverlaps(hybrids_tab,ref,by.x=c("chr.g","start.g","end.g"),by.y=c("seqnames","start","end"),nomatch = 0)
    
    #check the overlap length of annotated genes
    overlapped[,overlap:=apply(overlapped[,.(start, end, start.g, end.g)],1,function(x) min(x[2],x[4])-max(x[1],x[3])+1)]
    
    res <- overlapped[, .(paste0(unique(.SD[overlap == max(overlap),gene_name]),collapse = ","),paste0(unique(.SD[overlap == max(overlap),gene_biotype]),collapse = ",")), by=.(chr.g,start.g,end.g)]
    #resolve targets with overlapping annotation
    #add read name
    hybrids_tab <- merge(res, hybrids_tab, by = c("chr.g","start.g","end.g"), all.y = T)
    hybrids_tab[,gene_name := V1]
    hybrids_tab[,gene_biotype := V2]
    hybrids_tab$V1 <- NULL
    hybrids_tab$V2 <- NULL
    
    #add annotation of genomic feature location
    #########################
    annot_feat <- as.data.table(annot_feat)
    annot_feat <- annot_feat[,.(seqnames, start, end, type)]
    names(annot_feat) <- c("seqnames","start","end","loc")
    setkeyv(annot_feat, c("seqnames","start","end"))
    #add end position of target alignment
    setkeyv(hybrids_tab, c("chr.g","start.g","end.g"))
    
    
    overlapped_feat2 <- foverlaps(hybrids_tab, annot_feat, by.x=c("chr.g","start.g","end.g"), by.y=c("seqnames","start","end"))
    
    res_annot <- overlapped_feat2[, .(paste0(unique(loc),collapse = ",")), by=.(chr.g,start.g,end.g)]
    res_annot2 <- merge(res_annot, overlapped_feat2,by = c("chr.g","start.g","end.g"), all.y = T)
    
    res_annot2[, feature := V1]
    res_annot2$start <- NULL
    res_annot2$end <- NULL
    res_annot2$V1 <- NULL
    res_annot2$loc <- NULL
    hybrids_tab <- unique(res_annot2)
    #########################################
    #fix mirna names by real sequence
    
    fixed_smallRNA_names <- hybrids_tab[,paste0(unique(unlist(strsplit(noncodingRNA,"\\|"))),collapse = "|"),by = .(noncodingRNA_real_seq)]
    hybrids_tab <- merge(hybrids_tab,fixed_smallRNA_names, by ="noncodingRNA_real_seq")
    hybrids_tab[, noncodingRNA := V1]
    hybrids_tab$V1 <- NULL
    if (is_umi == TRUE){
      hybrids_tab <- hybrids_tab[,.(chr.g,start.g,end.g,gene_name,gene_biotype,noncodingRNA,noncodingRNA_fam,forward_dir,Ndups,Nunique,flag.g,strand.g,noncodingRNA_nomm,mrna_length,cov,
                                    Nreads,mrna_aln,mrna_aln2,mrna_aln5,mrna_aln10,mrna_aln20,noncodingRNA_seq,noncodingRNA_type,repeatmasker,repeat_fam,feature,noncodingRNA_real_seq,
                                    rep,dups_perc,seq.g)]
    } else {
      hybrids_tab <- hybrids_tab[,.(chr.g,start.g,end.g,gene_name,gene_biotype,noncodingRNA,noncodingRNA_fam,forward_dir,Ndups,flag.g,strand.g,noncodingRNA_nomm,mrna_length,cov,
                                    Nreads,mrna_aln,mrna_aln2,mrna_aln5,mrna_aln10,mrna_aln20,noncodingRNA_seq,noncodingRNA_type,repeatmasker,repeat_fam,feature,noncodingRNA_real_seq,
                                    rep,seq.g)]
    }

    # write_cofold_merged
    all_reps_collap <- copy(hybrids_tab)
    all_reps_collap_sort <- all_reps_collap[order(-Nunique)]
    seqs <- sapply(seq_along(all_reps_collap_sort$noncodingRNA), function(x) paste0(all_reps_collap_sort[x,seq.g],"&",if(all_reps_collap_sort[x,noncodingRNA_type == "miRNA" & !is.na(noncodingRNA_seq)]){all_reps_collap_sort[x,noncodingRNA_seq]}else{all_reps_collap_sort[x,noncodingRNA_real_seq]}))
    fa <- character(2 * length(seqs))
    fa[c(TRUE, FALSE)] = sprintf(">%s", seqs)
    fa[c(FALSE, TRUE)] = seqs
    writeLines(fa,paste0(c("hyb_pairs/cofold_pairs_collapsed_merged_",is_ambiguous,"/Merged_all_replicates_filtered_collapsed_hybrid_sequence_",type,".fa"),collapse = ""))
    system(paste0(c("( cd hyb_pairs/cofold_pairs_collapsed_merged_",is_ambiguous,"/all_types/", " && RNAcofold --noLP --output-format='D' < ../Merged_all_replicates_filtered_collapsed_hybrid_sequence_",type,".fa > ../Merged_",type,".collapsed_cofold.csv )"),collapse = ""))
    
    # add_cofold_and_unified 
    cofold <- fread(paste0(c("hyb_pairs/cofold_pairs_collapsed_merged_",is_ambiguous,"/Merged_",type,".collapsed_cofold.csv"),collapse = ""))
    cofold[,noncodingRNA := tstrsplit(seq_id,"&")[2]]
    if(type == "mirna"){
      all_reps_collap_sort[, A_perc := round(stri_count(noncodingRNA_seq, fixed = "A")/nchar(noncodingRNA_seq)*100,1)]
      all_reps_collap_sort[, C_perc := round(stri_count(noncodingRNA_seq, fixed = "C")/nchar(noncodingRNA_seq)*100,1)]
      all_reps_collap_sort[, G_perc := round(stri_count(noncodingRNA_seq, fixed = "G")/nchar(noncodingRNA_seq)*100,1)]
      all_reps_collap_sort[, T_perc := round(stri_count(noncodingRNA_seq, fixed = "T")/nchar(noncodingRNA_seq)*100,1)]
    }else{
      all_reps_collap_sort[, A_perc := round(stri_count(noncodingRNA_real_seq, fixed = "A")/nchar(noncodingRNA_real_seq)*100,1)]
      all_reps_collap_sort[, C_perc := round(stri_count(noncodingRNA_real_seq, fixed = "C")/nchar(noncodingRNA_real_seq)*100,1)]
      all_reps_collap_sort[, G_perc := round(stri_count(noncodingRNA_real_seq, fixed = "G")/nchar(noncodingRNA_real_seq)*100,1)]
      all_reps_collap_sort[, T_perc := round(stri_count(noncodingRNA_real_seq, fixed = "T")/nchar(noncodingRNA_real_seq)*100,1)]
    }
    
    all_reps_collap_sort<- cbind(all_reps_collap_sort,cofold[,.(mfe_struct,mfe)])
    # convert to the energy numbers according their formula
    all_reps_collap_sort[,ene :=5*(mfe-(-11))/(-5) ]
    all_reps_collap_sort[,ene := fifelse(ene >= 5,5,ene)]
    all_reps_collap_sort[,ene := fifelse(ene <= 0,0,ene)]
    if (is_umi == TRUE) {
      all_reps_collap_sort <- all_reps_collap_sort[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, Nunique,
                                                       dups_perc, repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                       mfe, ene, noncodingRNA_seq,rep)] 
      all_reps_collap_sort <- all_reps_collap_sort[order(-Nunique)]
    } else {
      all_reps_collap_sort <- all_reps_collap_sort[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, 
                                                       repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                       mfe, ene, noncodingRNA_seq,rep)] 
      all_reps_collap_sort <- all_reps_collap_sort[order(-Ndups)]
    }
    
    #resolve replicates order
    all_reps_collap_sort[, rep := sapply(seq_along(rep), function(x) {paste(sort(unique(unlist(tstrsplit(rep[x],",")))),collapse = ",")})]
    all_reps_collap_sort[, noncodingRNA := sapply(seq_along(noncodingRNA), function(x) paste0(unique(tstrsplit(noncodingRNA[x],"\\|")),collapse = "|"))]
    write.table(all_reps_collap_sort,paste0(c("hyb_pairs/Merged.hybrids_deduplicated_filtered_collapsed_",type,"_",is_ambiguous,".tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    return(all_reps_collap_sort)
  }
  
  if (is_ambiguous == "ambiguous"){
    all_reps_tab_collapsed <- merge_and_collapse(all_reps_tab,repeatmasker,ref,annot_feat,"all")
  } else {
    for (type in names(table(all_reps_tab$noncodingRNA_type))){
      
      all_reps_collapsed <- merge_and_collapse(all_reps_tab,repeatmasker,ref,annot_feat,type)
      if (!exists("all_reps_tab_collapsed")){
        all_reps_tab_collapsed <- copy(all_reps_collapsed)
      } else {
        all_reps_tab_collapsed <- rbind(all_reps_tab_collapsed,all_reps_collapsed)
      }
      
    }
  }
  write.table(all_reps_tab_collapsed,paste0(c("hyb_pairs/Merged.hybrids_deduplicated_filtered_collapsed_",is_ambiguous,".tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
}

# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)  