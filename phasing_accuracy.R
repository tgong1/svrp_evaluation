#!/usr/bin/env Rscript
#

bedpe_forbedtools <- function(input, bkpt_diff_threshold){
  options(scipen=999)
  diff <- bkpt_diff_threshold/2
  input$SVLEN <- as.numeric(input$SVLEN)
  input$END <- as.numeric(input$END)
  bedpe <- data.frame(chrom1 = input$chrom,
                      start1 = ifelse(input$POS-diff <0, 0,input$POS-diff-1),
                      end1 = input$POS+diff,
                      chrom2 = input$chrom,
                      start2 = ifelse(input$END-diff <0, 0, input$END-diff-1),
                      end2 = input$END+diff,
                      ID = input$ID,
                      SVTYPE=input$SVTYPE)
  return(bedpe)
}

compare_truth_bedtools <- function(eval, truth, DIR, bkpt_diff_threshold, out_name){
  bedpe <- bedpe_forbedtools(eval, bkpt_diff_threshold)
  write.table(bedpe, paste0(DIR,"eval.sv.bedpe"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  bedpe <- bedpe_forbedtools(truth, bkpt_diff_threshold)
  write.table(bedpe, paste0(DIR,"truth.sv.bedpe"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  
  system(paste0("bedtools pairtopair -a ",DIR,"eval.sv.bedpe"," -b ",DIR,"truth.sv.bedpe > ",DIR, out_name ))
  if(file.info(paste0(DIR, out_name))$size==0){
    overlap <- c()
  }else{
    overlap <- read.table(paste0(DIR, out_name), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    colnames(overlap) <- c("SV_chrom1", "SV_start1","SV_end1",
                           "SV_chrom2", "SV_start2","SV_end2",
                           "SV_name","SVTYPE",
                           "hf_chrom1","hf_start1","hf_end1",
                           "hf_chrom2","hf_start2","hf_end2",
                           "hf_name","hf_SVTYPE",
                           "overlap1", "overlap2")
    overlap$diff1 <- abs(overlap$hf_start1- overlap$SV_start1)
    overlap$diff2 <- abs(overlap$hf_start2- overlap$SV_start2)
    overlap$diff <- overlap$diff1 + overlap$diff2 
    
    overlap <- overlap[overlap$SVTYPE==overlap$hf_SVTYPE,]
    overlap <- overlap[!duplicated(overlap),]
  }
  
  return(overlap)
}

remove_duplicated_overlap <- function(overlap){
  duplicated <- overlap[duplicated(overlap$SV_name),]
  if(nrow(duplicated)!=0){
    unique <- c()
    for(i in c(1: nrow(duplicated))){
      tmp <- overlap[overlap$SV_name == duplicated[i,]$SV_name,]
      unique <- rbind(unique, tmp[which(tmp$diff ==min(tmp$diff)),])
    }
    
    overlap <- rbind(overlap[!(overlap$SV_name %in% duplicated$SV_name),],
                     unique)
  }
  
  return(overlap)
}


generate_phased_SNV_SV <- function(eval_SV, eval_SNV, truth_SV, DIR, bkpt_diff_threshold, out_name){
  overlap <- compare_truth_bedtools(eval = eval_SV, truth = truth_SV, DIR, bkpt_diff_threshold, out_name)
  if(is.null(overlap)){
    eval_SV$new_ID <- eval_SV$ID
  }else{
    overlap_unique <- remove_duplicated_overlap(overlap)
    eval_SV$new_ID <- eval_SV$ID
    vs_index <- match(eval_SV$ID, overlap_unique$SV_name)
    index <- vs_index[!is.na(vs_index)]
    eval_SV$new_ID[!is.na(vs_index)] <- overlap_unique[index,]$hf_name
  }
  longphased_phased_SNV_SV <- rbind(eval_SNV, eval_SV)
  df <- longphased_phased_SNV_SV
  chrOrder <- paste0("chr",c(1:22))
  df$chrom <- factor(df$chrom, chrOrder, ordered=TRUE)
  df <- df[do.call(order, df[, c( "chrom", "POS")]), ]
  
  return(df)
}

compare_eval_truth_GT <- function(eval_ALL, truth_ALL){
  eval_ALL$pair_index <- match(eval_ALL$new_ID, truth_ALL$ID)
  eval_ALL$truth_ID <- eval_ALL$pair_index
  index <- eval_ALL$pair_index[!is.na(eval_ALL$pair_index)]
  eval_ALL$truth_ID[!is.na(eval_ALL$pair_index)] <- truth_ALL[index,]$ID
  eval_ALL$truth_GT[!is.na(eval_ALL$pair_index)] <- truth_ALL[index,]$GT
  
  longphased_phased_match_truth <- eval_ALL[(!is.na(eval_ALL$pair_index)),]
  phased_het_SV_index <- which(grepl("SV", longphased_phased_match_truth$ID))
  longphased_phased_match_truth_SV <- longphased_phased_match_truth[phased_het_SV_index,]
  if(length(phased_het_SV_index) != 0){
    count <- c()
    POS_dist <- c()
    for(i in c(1:length(phased_het_SV_index))){
      tmp <- longphased_phased_match_truth[c(phased_het_SV_index[i]-1, phased_het_SV_index[i]),]
      tmp_count <- sum(sum(tmp$GT == tmp$truth_GT)==1)
      count <- c(count, tmp_count)
      tmp_dist <- abs(tmp$POS[1]- tmp$POS[2])
      POS_dist <- c(POS_dist, tmp_dist)
    }
    longphased_phased_match_truth_SV$count_switch1 <- count
    longphased_phased_match_truth_SV$POS_distance1 <- POS_dist
    
    count <- c()
    POS_dist <- c()
    for(i in c(1:length(phased_het_SV_index))){
      tmp <- longphased_phased_match_truth[c(phased_het_SV_index[i], phased_het_SV_index[i]+1),]
      tmp_count <- sum(sum(tmp$GT == tmp$truth_GT)==1)
      count <- c(count, tmp_count)
      
      tmp_dist <- abs(tmp$POS[1]- tmp$POS[2])
      POS_dist <- c(POS_dist, tmp_dist)
    }
    longphased_phased_match_truth_SV$count_switch2 <- count
    longphased_phased_match_truth_SV$POS_distance2 <- POS_dist
  }
  
  return(longphased_phased_match_truth_SV)
}

flip_rate_each_chrom_LCL6 <- function(eval_SV,eval_SNV, truth_SV, threshold, i, DIR, DIR_truth){
  longphased_phased_SNV <- read.table(paste0(DIR, eval_SNV), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(longphased_phased_SNV) <- c("chrom","POS","ID","REF","ALT","FILTER","GT")
  longphased_phased_SNV <- cbind(longphased_phased_SNV[,c(1:6)],
                                 SVTYPE=".",
                                 SVLEN=".",
                                 END=".",
                                 GT=longphased_phased_SNV[,7],
                                 new_ID =longphased_phased_SNV$ID)
  ###
  truth <- paste0(DIR_truth,"LCL6.high.confidence_SNV_SV_pedigree_phased_IDmod_chr",i,"_sorted.txt")
  truth_ALL <- read.table(truth, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(truth_ALL) <- c("chrom","POS","ID","REF","ALT","FILTER","GT")
  
  ###
  longphased_phased_SNV_SV_sorted <- generate_phased_SNV_SV(eval_SV = eval_SV, 
                                                            eval_SNV = longphased_phased_SNV,
                                                            truth_SV = truth_SV, 
                                                            DIR, bkpt_diff_threshold = threshold, 
                                                            out_name = "LCL6longphase_vs_hf.bedpe")
  longphased_phased_match_truth_SV <- compare_eval_truth_GT(eval_ALL = longphased_phased_SNV_SV_sorted,
                                                            truth_ALL)
  
  longphased_phased_match_truth_SV <- longphased_phased_match_truth_SV[(!is.na(longphased_phased_match_truth_SV$count_switch1)) &
                                                                         (!is.na(longphased_phased_match_truth_SV$count_switch2)), ]
  df <- c()
  for(SVTYPE in c("DEL","INS")){
    result <- longphased_phased_match_truth_SV[longphased_phased_match_truth_SV$SVTYPE==SVTYPE,]
    N_flip <- sum(result$count_switch1==1 &
                    result$count_switch2==1)
    N_zero_switch <- sum(result$count_switch1==0 &
                           result$count_switch2==0)
    N_one_switch <- sum((result$count_switch1==1 &
                           result$count_switch2==0)|
                          result$count_switch1==0 &
                          result$count_switch2==1)
    
    df <- data.frame(rbind(df, c(chrom=paste0("chr",i),
                                 paired_SV_count = nrow(result),
                                 N_flip = N_flip,
                                 N_zero_switch = N_zero_switch,
                                 N_one_switch = N_one_switch,
                                 SVTYPE = SVTYPE)))
    
  }
  
  return(df)
}

read_LCL6_longphased_phased_SNV_withBlock <- function(eval_SNV, PS_ID, DIR){
  longphased_phased_SNV <- read.table(paste0(DIR, eval_SNV), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(longphased_phased_SNV) <- c("chrom","POS","ID","REF","ALT","FILTER","GT", "PS")
  longphased_phased_SNV <- longphased_phased_SNV[longphased_phased_SNV$PS == PS_ID,]
  
  if(nrow(longphased_phased_SNV) != 0){
    longphased_phased_SNV <- cbind(longphased_phased_SNV[,c(1:6)],
                                   SVTYPE=".",
                                   SVLEN=".",
                                   END=".",
                                   GT=longphased_phased_SNV[,7],
                                   new_ID =longphased_phased_SNV$ID)}
  
  return(longphased_phased_SNV)
}

flip_rate_each_chrom_LCL6_withBlock <- function(eval_SV, longphased_phased_SNV, truth_SV, threshold, i, DIR_truth){
  ###
  truth <- paste0(DIR_truth, "LCL6.high.confidence_SNV_SV_pedigree_phased_IDmod_chr",i,"_sorted.txt")
  truth_ALL <- read.table(truth, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(truth_ALL) <- c("chrom","POS","ID","REF","ALT","FILTER","GT")
  
  ###
  longphased_phased_SNV_SV_sorted <- generate_phased_SNV_SV(eval_SV = eval_SV, 
                                                            eval_SNV = longphased_phased_SNV,
                                                            truth_SV = truth_SV, 
                                                            DIR, bkpt_diff_threshold = threshold, 
                                                            out_name = "LCL6longphase_vs_hf.bedpe")
  longphased_phased_match_truth_SV <- compare_eval_truth_GT(eval_ALL = longphased_phased_SNV_SV_sorted,
                                                            truth_ALL)
  
  longphased_phased_match_truth_SV <- longphased_phased_match_truth_SV[(!is.na(longphased_phased_match_truth_SV$count_switch1)) &
                                                                         (!is.na(longphased_phased_match_truth_SV$count_switch2)), ]
  df <- c()
  for(SVTYPE in c("DEL","INS")){
    result <- longphased_phased_match_truth_SV[longphased_phased_match_truth_SV$SVTYPE==SVTYPE,]
    N_flip <- sum(result$count_switch1==1 &
                    result$count_switch2==1)
    N_zero_switch <- sum(result$count_switch1==0 &
                           result$count_switch2==0)
    N_one_switch <- sum((result$count_switch1==1 &
                           result$count_switch2==0)|
                          result$count_switch1==0 &
                          result$count_switch2==1)
    
    df <- data.frame(rbind(df, c(chrom=paste0("chr",i),
                                 paired_SV_count = nrow(result),
                                 N_flip = N_flip,
                                 N_zero_switch = N_zero_switch,
                                 N_one_switch = N_one_switch,
                                 SVTYPE = SVTYPE)))
    
  }
  
  return(df)
}

SVTYPE_flip_rate_HG002 <- function(longphased_phased_match_truth_SV, threshold,df_sum){
  
  for(SVTYPE in c("DEL","INS")){
    result <- longphased_phased_match_truth_SV[longphased_phased_match_truth_SV$SVTYPE==SVTYPE,]
    N_flip <- sum(result$count_switch1==1 &
                    result$count_switch2==1)
    N_zero_switch <- sum(result$count_switch1==0 &
                           result$count_switch2==0)
    N_one_switch <- sum(result$count_switch1==0 & result$count_switch2==1)+ sum(result$count_switch1==1 & result$count_switch2==0)
    #Compute an SV flip rate by dividing the number of SVs with 1 flip by the total number of assessed HET SVs (total no. of assessed SVs = (SVs with 1 flip + SVs with no flips and no switch errors)) in a given sample. 
    
    tmp <- data.frame(N_total = nrow(result),
                      N_flip = N_flip, 
                      N_one_switch = N_one_switch, 
                      N_zero_switch = N_zero_switch, 
                      flip_rate = N_flip/(N_flip+ N_zero_switch),
                      threshold = threshold,
                      SVTYPE= SVTYPE)
    df_sum <- rbind(df_sum, tmp)
  }
  return(df_sum)
}

SVTYPE_flip_rate_HG002_withBlock <- function(longphased_phased_match_truth_SV, threshold){
  df_sum <- c()
  for(SVTYPE in c("DEL","INS")){
    result <- longphased_phased_match_truth_SV[longphased_phased_match_truth_SV$SVTYPE==SVTYPE,]
    N_flip <- sum(result$count_switch1==1 &
                    result$count_switch2==1)
    N_zero_switch <- sum(result$count_switch1==0 &
                           result$count_switch2==0)
    N_one_switch <- sum(result$count_switch1==0 & result$count_switch2==1)+ sum(result$count_switch1==1 & result$count_switch2==0)
    #Compute an SV flip rate by dividing the number of SVs with 1 flip by the total number of assessed HET SVs (total no. of assessed SVs = (SVs with 1 flip + SVs with no flips and no switch errors)) in a given sample. 
    
    tmp <- data.frame(N_total = nrow(result),
                      N_flip = N_flip, 
                      N_one_switch = N_one_switch, 
                      N_zero_switch = N_zero_switch, 
                      flip_rate = N_flip/(N_flip+ N_zero_switch),
                      threshold = threshold,
                      SVTYPE= SVTYPE)
    df_sum <- rbind(df_sum, tmp)
  }
  return(df_sum)
}

flip_rate_each_chrom_HG002 <- function(eval_SV, truth_ALL,truth_SV, eval_SNV_file, threshold, DIR){
  shapeit5_phased_SNV <- read.table(paste0(DIR, eval_SNV_file), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(shapeit5_phased_SNV) <- c("chrom","POS","ID","REF","ALT","FILTER","GT")
  shapeit5_phased_SNV <- cbind(shapeit5_phased_SNV[,c(1:6)],
                               SVTYPE=".",
                               SVLEN=".",
                               END=".",
                               GT=shapeit5_phased_SNV[,7],
                               new_ID =shapeit5_phased_SNV$ID)
  
  ###
  longphased_phased_SNV_SV_sorted <- generate_phased_SNV_SV(eval_SV = eval_SV[eval_SV$chrom==paste0("chr",i),],
                                                            eval_SNV = shapeit5_phased_SNV,
                                                            truth_SV = truth_SV[truth_SV$chrom==paste0("chr",i),], 
                                                            DIR, bkpt_diff_threshold=threshold, 
                                                            out_name = "HG002shapeit5_vs_hf.bedpe")
  
  longphased_phased_match_truth_SV <- compare_eval_truth_GT(eval_ALL = longphased_phased_SNV_SV_sorted,
                                                            truth_ALL = truth_ALL[truth_ALL$chrom==paste0("chr",i),])
  
  longphased_phased_match_truth_SV <- longphased_phased_match_truth_SV[(!is.na(longphased_phased_match_truth_SV$count_switch1)) &
                                                                         (!is.na(longphased_phased_match_truth_SV$count_switch2)), ]
  df <- c()
  for(SVTYPE in c("DEL","INS")){
    result <- longphased_phased_match_truth_SV[longphased_phased_match_truth_SV$SVTYPE==SVTYPE,]
    N_flip <- sum(result$count_switch1==1 &
                    result$count_switch2==1)
    N_zero_switch <- sum(result$count_switch1==0 &
                           result$count_switch2==0)
    N_one_switch <- sum((result$count_switch1==1 &
                           result$count_switch2==0)|
                          result$count_switch1==0 &
                          result$count_switch2==1)
    
    df <- data.frame(rbind(df, c(chrom=paste0("chr",i),
                                 paired_SV_count = nrow(result),
                                 N_flip = N_flip,
                                 N_zero_switch = N_zero_switch,
                                 N_one_switch = N_one_switch,
                                 SVTYPE = SVTYPE)))
    
  }
  
  return(df)
}