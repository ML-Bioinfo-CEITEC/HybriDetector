library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

run_all <- function(args){
  hybrids <- args[1]
  sample <- args[2]
  repeatmasker_bed <-args[3]
  transcripts <- args[4]
  features <- args[5]
  is_umi <- args[6]
  merged <- args[7]
  is_ambiguous <- args[8]

  input_hybrids <- fread(hybrids)
  repeatmasker <- fread(repeatmasker_bed)
  ref <- as.data.table(rtracklayer::import(transcripts))[type == "transcript", c("seqnames","start","end","gene_name","strand","gene_biotype"), with=F]
  annot_feat <- rtracklayer::import(features)

  input_hybrids[, new_mid := round((start.g+end.g)/2,0)]
  input_hybrids[, new_start := new_mid -24]
  input_hybrids[, new_end := new_mid + 25]
  input_hybrids_new_pos <- input_hybrids[,c(1,(dim(input_hybrids)[2]-1),(dim(input_hybrids)[2]),4:(dim(input_hybrids)[2]-3)),with = F] 
  
  names(input_hybrids_new_pos)[names(input_hybrids_new_pos) == 'new_start'] <- 'start.g'
  names(input_hybrids_new_pos)[names(input_hybrids_new_pos) == 'new_end'] <- 'end.g'
  
  chromosome <- paste0("chr",input_hybrids_new_pos$chr.g)
  chromosome <- gsub("chrMT","chrM",chromosome)
  
  ranges <- GRanges( chromosome,
                     IRanges(start=input_hybrids_new_pos$start.g, end = input_hybrids_new_pos$end.g),
                     strand=input_hybrids_new_pos$strand.g )
  target_sequences <- as.data.table(Views(Hsapiens,ranges))
  input_hybrids_new_pos[,seq.g := target_sequences$dna]


  seqs <- sapply(seq_along(input_hybrids_new_pos$noncodingRNA), function(x) paste0(input_hybrids_new_pos[x,seq.g],"&",if(input_hybrids_new_pos[x,noncodingRNA_type == "miRNA" & !is.na(noncodingRNA_seq)]){input_hybrids_new_pos[x,noncodingRNA_seq]}else{input_hybrids_new_pos[x,noncodingRNA_real_seq]}))
  fa <- character(2 * length(seqs))
  fa[c(TRUE, FALSE)] = sprintf(">%s", seqs)
  fa[c(FALSE, TRUE)] = seqs
  writeLines(fa, paste0(c("hyb_pairs/cofold_pairs_collapsed_unified_",is_ambiguous,"/",sample,".filtered_collapsed_unified_hybrid_sequence.fa"),collapse = ""))
  system(paste0(c("( cd hyb_pairs/cofold_pairs_collapsed_unified_",is_ambiguous,"/",sample,"/", " && RNAcofold --noLP --output-format='D' < ../",sample,".filtered_collapsed_unified_hybrid_sequence.fa > ../",sample,".collapsed_unified_cofold.csv )"),collapse = ""))
  
  cofold <- fread(paste0(c("hyb_pairs/cofold_pairs_collapsed_unified_",is_ambiguous,"/",sample,".collapsed_unified_cofold.csv"),collapse = ""))
  cofold[,rna := tstrsplit(seq_id,"&")[1]]
  cofold[,noncodingRNA := tstrsplit(seq_id,"&")[2]]
  input_hybrids_new_pos$mfe_struct <- NULL
  input_hybrids_new_pos$mfe <- NULL
  input_hybrids_new_pos$ene <- NULL
  input_hybrids_new_pos <- cbind(input_hybrids_new_pos,cofold[,c(4,5)])
  input_hybrids_new_pos[,ene :=5*(mfe-(-11))/(-5) ]
  input_hybrids_new_pos[,ene := fifelse(ene >= 5,5,ene)]
  input_hybrids_new_pos[,ene := fifelse(ene <= 0,0,ene)]
  input_hybrids_new_pos <- input_hybrids_new_pos[,c(1:(dim(input_hybrids_new_pos)[2]-5),(dim(input_hybrids_new_pos)[2]-2):dim(input_hybrids_new_pos)[2],(dim(input_hybrids_new_pos)[2]-4),(dim(input_hybrids_new_pos)[2]-3)),with = F]
  if (is_umi == TRUE) {
    input_hybrids_new_pos <- input_hybrids_new_pos[order(-Nunique)]
  } else {
    input_hybrids_new_pos <- input_hybrids_new_pos[order(-Ndups)]
    
  }
  #####################################
  #check repeat_masker
  #add end position for repeatmasker overlap
  names(repeatmasker) <- c("seqnames","start","end","repeatmasker","fam","type")
  repeatmasker <- repeatmasker[,.(seqnames,start,end,repeatmasker,fam)]
  setkeyv(repeatmasker, c("seqnames","start","end"))
  setkeyv(input_hybrids_new_pos, c("chr.g","start.g","end.g"))
  
  overlapped_repeatmasker <- foverlaps(input_hybrids_new_pos, repeatmasker, by.x=c("chr.g","start.g","end.g"), by.y=c("seqnames","start","end"))
  
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
  input_hybrids_new_pos <- unique(res_annot2)

  #########################################
  #fix mirna names by real sequence

  fixed_smallRNA_names <- input_hybrids_new_pos[,paste0(unique(unlist(strsplit(noncodingRNA,"\\|"))),collapse = "|"),by = .(noncodingRNA_real_seq)]
  input_hybrids_new_pos <- merge(input_hybrids_new_pos,fixed_smallRNA_names, by ="noncodingRNA_real_seq")
  input_hybrids_new_pos[, noncodingRNA := V1]
  input_hybrids_new_pos$V1 <- NULL
  names(input_hybrids_new_pos)[names(input_hybrids_new_pos) == 'smallRNA_fam'] <- 'noncodingRNA_fam'
  names(input_hybrids_new_pos)[names(input_hybrids_new_pos) == 'smallRNA_nomm'] <- 'noncodingRNA_nomm'
  
  if (is_umi == TRUE) {
    if (merged == TRUE) {
      input_hybrids_new_pos <- input_hybrids_new_pos[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, Nunique,
                                                         dups_perc, repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                         mfe, ene, noncodingRNA_seq,rep)] 
      all_hybrids_collapsed_unified_pub <- input_hybrids_new_pos[order(-Nunique)]
      all_hybrids_collapsed_unified_pub[, ID := seq_along(gene_name)]
      all_hybrids_collapsed_unified_pub <- all_hybrids_collapsed_unified_pub[,.(ID, seq.g, noncodingRNA_real_seq, forward_dir, Ndups, Nunique, rep, Nreads, cov, chr.g, start.g, end.g, strand.g, gene_name, gene_biotype, feature, repeatmasker, repeat_fam, noncodingRNA, 
                                           noncodingRNA_type, noncodingRNA_fam,noncodingRNA_nomm,mfe_struct,mfe)]
      names(all_hybrids_collapsed_unified_pub) <- c("ID", "Genomic fragment sequence", "Driver fragment sequence", "Target first in chimera", "N chimeric reads", "N deduplicated chimeric reads", "Chimera in replicates", "N non-chimeric reads", "Coverage non-chimeric reads", "Chromosome", "Start", "End", "Strand", "Overlapping gene name", "Overlapping gene biotype", "Overlapping gene feature", "Overlapping repeatmasker", "Overlapping repeatmasker family", "Driver name", "Driver type", "miRNA Family", "Driver alignment no mismatch", "Cofold structure", "Cofold MFE")
      
    } else {
      input_hybrids_new_pos <- input_hybrids_new_pos[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, Nunique,
                                                         dups_perc, repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                          mfe, ene, noncodingRNA_seq)] 
      all_hybrids_collapsed_unified_pub <- input_hybrids_new_pos[order(-Nunique)]
      all_hybrids_collapsed_unified_pub[, ID := seq_along(gene_name)]
      all_hybrids_collapsed_unified_pub <- all_hybrids_collapsed_unified_pub[,.(ID, seq.g, noncodingRNA_real_seq, forward_dir, Ndups, Nunique, Nreads, cov, chr.g, start.g, end.g, strand.g, gene_name, gene_biotype, feature, repeatmasker, repeat_fam, noncodingRNA, 
                                                                                noncodingRNA_type, noncodingRNA_fam,noncodingRNA_nomm,mfe_struct,mfe)]
      names(all_hybrids_collapsed_unified_pub) <- c("ID", "Genomic fragment sequence", "Driver fragment sequence", "Target first in chimera", "N chimeric reads", "N deduplicated chimeric reads", "N non-chimeric reads", "Coverage non-chimeric reads", "Chromosome", "Start", "End", "Strand", "Overlapping gene name",  "Overlapping gene biotype", "Overlapping gene feature", "Overlapping repeatmasker", "Overlapping repeatmasker family", "Driver name", "Driver type", "miRNA Family", "Driver alignment no mismatch", "Cofold structure", "Cofold MFE")
      
    }
    all_hybrids_collapsed_unified <- input_hybrids_new_pos[order(-Nunique)]
        
  } else {
    if (merged == TRUE) {
      input_hybrids_new_pos <- input_hybrids_new_pos[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, 
                                                         repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                         mfe, ene, noncodingRNA_seq,rep)]
      all_hybrids_collapsed_unified_pub <- input_hybrids_new_pos[order(-Ndups)]
      all_hybrids_collapsed_unified_pub[, ID := seq_along(gene_name)]
      all_hybrids_collapsed_unified_pub <- all_hybrids_collapsed_unified_pub[,.(ID, seq.g, noncodingRNA_real_seq, forward_dir, Ndups, rep, Nreads, cov, chr.g, start.g, end.g, strand.g, gene_name, gene_biotype, feature, repeatmasker, repeat_fam, noncodingRNA, 
                                                                                noncodingRNA_type, noncodingRNA_fam,noncodingRNA_nomm,mfe_struct,mfe)]
      names(all_hybrids_collapsed_unified_pub) <- c("ID", "Genomic fragment sequence", "Driver fragment sequence", "Target first in chimera", "N chimeric reads", "Chimera in replicates", "N non-chimeric reads", "Coverage non-chimeric reads", "Chromosome", "Start", "End", "Strand", "Overlapping gene name",  "Overlapping gene biotype", "Overlapping gene feature", "Overlapping repeatmasker", "Overlapping repeatmasker family", "Driver name", "Driver type", "miRNA Family", "Driver alignment no mismatch", "Cofold structure", "Cofold MFE")
      
    } else {
      input_hybrids_new_pos <- input_hybrids_new_pos[ ,.(chr.g,start.g,end.g,flag.g,strand.g,gene_name,gene_biotype,feature,seq.g,noncodingRNA_real_seq,noncodingRNA, noncodingRNA_type, noncodingRNA_fam, forward_dir, Ndups, 
                                                         repeatmasker, repeat_fam, noncodingRNA_nomm, mrna_length, cov, Nreads, mrna_aln, mrna_aln2, mrna_aln5, mrna_aln10, mrna_aln20, A_perc, C_perc, G_perc, T_perc, mfe_struct,
                                                         mfe, ene, noncodingRNA_seq)] 
      all_hybrids_collapsed_unified_pub <- input_hybrids_new_pos[order(-Ndups)]
      all_hybrids_collapsed_unified_pub[, ID := seq_along(gene_name)]
      all_hybrids_collapsed_unified_pub <- all_hybrids_collapsed_unified_pub[,.(ID, seq.g, noncodingRNA_real_seq, forward_dir, Ndups, Nreads, cov, chr.g, start.g, end.g, strand.g, gene_name, gene_biotype, feature, repeatmasker, repeat_fam, noncodingRNA, 
                                                                                noncodingRNA_type, noncodingRNA_fam,noncodingRNA_nomm,mfe_struct,mfe)]
      names(all_hybrids_collapsed_unified_pub) <- c("ID", "Genomic fragment sequence", "Driver fragment sequence", "Target first in chimera", "N chimeric reads", "N non-chimeric reads", "Coverage non-chimeric reads", "Chromosome", "Start", "End", "Strand", "Overlapping gene name", "Overlapping gene biotype", "Overlapping gene feature", "Overlapping repeatmasker", "Overlapping repeatmasker family", "Driver name", "Driver type", "miRNA Family", "Driver alignment no mismatch", "Cofold structure", "Cofold MFE")
      
    }
    all_hybrids_collapsed_unified <- input_hybrids_new_pos[order(-Ndups)]
  }

  #############
  all_hybrids_collapsed_unified_high_confidence <- all_hybrids_collapsed_unified[noncodingRNA_nomm == TRUE]
  all_hybrids_collapsed_unified_pub_high_confidence <- all_hybrids_collapsed_unified_pub[`Driver alignment no mismatch` == TRUE]
  if (merged == TRUE) {
    write.table(all_hybrids_collapsed_unified,paste0(c("hyb_pairs/Merged.unified_length_all_types_",is_ambiguous,".tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_high_confidence,paste0(c("hyb_pairs/Merged.unified_length_all_types_",is_ambiguous,"_high_confidence_preclusters.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_pub,paste0(c("hyb_pairs/Merged.unified_length_all_type_",is_ambiguous,"_finalout.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_pub_high_confidence,paste0(c("hyb_pairs/Merged.unified_length_all_types_",is_ambiguous,"_high_confidence_finalout.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
  } else {
    write.table(all_hybrids_collapsed_unified,paste0(c("hyb_pairs/",sample,".unified_length_all_types_",is_ambiguous,".tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_high_confidence,paste0(c("hyb_pairs/",sample,".unified_length_all_types_",is_ambiguous,"_high_confidence_preclusters.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_pub,paste0(c("hyb_pairs/",sample,".unified_length_all_types_",is_ambiguous,"_finalout.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
    write.table(all_hybrids_collapsed_unified_pub_high_confidence,paste0(c("hyb_pairs/",sample,".unified_length_all_types_",is_ambiguous,"_high_confidence_pubout_finalout.tsv"),collapse = ""), sep = "\t",quote = F,row.names = F)
  }
}
# run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)


