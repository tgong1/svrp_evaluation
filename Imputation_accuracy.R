#!/usr/bin/env Rscript
#
imputation_TP_SNV_concordant_GT <- function(detected_SNV, GT_detected_SNV, imputed_SNV, GT_imputed_SNV, truth_ONTsamples){
  df_SNV <- c()
  for(i in c(1:nrow(truth_ONTsamples))){
    TP1 <- sum(detected_SNV[GT_detected_SNV[,i] %in% c("0/0","0|0"),]$ID %in% imputed_SNV[GT_imputed_SNV[,i]=="0|0",]$new_ID)
    TP2 <- sum(detected_SNV[GT_detected_SNV[,i] %in% c("0|1","1|0","0/1"),]$ID %in% imputed_SNV[GT_imputed_SNV[,i] %in% c("0|1","1|0"),]$new_ID)
    TP3 <- sum(detected_SNV[GT_detected_SNV[,i] %in% c("1/1","1|1"),]$ID %in% imputed_SNV[GT_imputed_SNV[,i]=="1|1",]$new_ID)
    n_total_detected <- nrow(detected_SNV[GT_detected_SNV[,i] %in% c("0/0","0|0"),]) + 
      nrow(detected_SNV[GT_detected_SNV[,i] %in% c("0|1","1|0","0/1"),]) + 
      nrow(detected_SNV[GT_detected_SNV[,i] %in% c("1/1","1|1"),])
    n_total_imputed <- nrow(imputed_SNV[GT_imputed_SNV[,i]=="0|0",]) + 
      nrow(imputed_SNV[GT_imputed_SNV[,i] %in% c("0|1","1|0"),]) + 
      nrow(imputed_SNV[GT_imputed_SNV[,i]=="1|1",])
    TP <- TP1+ TP2 + TP3
    recall <- TP/n_total_detected
    precision <- TP/n_total_imputed
    df_SNV <- rbind(df_SNV, data.frame(TP,
                                       truth= n_total_detected,
                                       imputed = n_total_imputed,
                                       recall, precision, 
                                       sample = truth_ONTsamples[i,1], 
                                       var_type ="SNV"))
  }
  return(df_SNV)
}

imputation_TP_indel_concordant_GT <- function(detected_indel, GT_detected_indel, imputed_indel, GT_imputed_indel, truth_ONTsamples){
  df_indel <- c()
  for(i in c(1:nrow(truth_ONTsamples))){
    TP1 <- sum(detected_indel[GT_detected_indel[,i] %in% c("0/0","0|0"),]$ID %in% imputed_indel[GT_imputed_indel[,i]=="0|0",]$new_ID)
    TP2 <- sum(detected_indel[GT_detected_indel[,i] %in% c("0|1","1|0","0/1"),]$ID %in% imputed_indel[GT_imputed_indel[,i] %in% c("0|1","1|0"),]$new_ID)
    TP3 <- sum(detected_indel[GT_detected_indel[,i] %in% c("1/1","1|1"),]$ID %in% imputed_indel[GT_imputed_indel[,i]=="1|1",]$new_ID)
    n_total_detected <- nrow(detected_indel[GT_detected_indel[,i] %in% c("0/0","0|0"),]) + 
      nrow(detected_indel[GT_detected_indel[,i] %in% c("0|1","1|0","0/1"),]) + 
      nrow(detected_indel[GT_detected_indel[,i] %in% c("1/1","1|1"),])
    n_total_imputed <- nrow(imputed_indel[GT_imputed_indel[,i]=="0|0",]) + 
      nrow(imputed_indel[GT_imputed_indel[,i] %in% c("0|1","1|0"),]) + 
      nrow(imputed_indel[GT_imputed_indel[,i]=="1|1",])
    
    TP <- TP1+ TP2 + TP3
    recall <- TP/n_total_detected
    precision <- TP/n_total_imputed
    df_indel <- rbind(df_indel, data.frame(TP, 
                                           truth= n_total_detected,
                                           imputed = n_total_imputed,
                                           recall, precision, 
                                           sample = truth_ONTsamples[i,1],
                                           var_type ="indel"))
  }
  return(df_indel)
}

imputation_TP_SV_concordant_GT <- function(detected_SV, GT_detected_SV, imputed_SV, GT_imputed_SV, truth_ONTsamples){
  df_SV <- c()
  for(i in c(1:nrow(truth_ONTsamples))){
    sample_truth1 = detected_SV[(GT_detected_SV[,i] %in% c("0/0","0|0")),]
    sample_imputed1 = imputed_SV[GT_imputed_SV[,i] == "0|0",]
    sample_TP1 <- sample_truth1[(sample_truth1$POS %in% as.vector(sapply(sample_imputed1$POS, function(x) (x-seq(-200,200))))) &
                                  sample_truth1$CHROM %in% sample_imputed1$CHROM &
                                  sample_truth1$SVTYPE %in% sample_imputed1$SVTYPE, ]
    
    sample_truth2 = detected_SV[(GT_detected_SV[,i] %in% c("0|1","1|0","0/1")),]
    sample_imputed2 = imputed_SV[GT_imputed_SV[,i] %in% c("0|1","1|0"),]
    sample_TP2 <- sample_truth2[(sample_truth2$POS %in% as.vector(sapply(sample_imputed2$POS, function(x) (x-seq(-200,200))))) &
                                  sample_truth2$CHROM %in% sample_imputed2$CHROM &
                                  sample_truth2$SVTYPE %in% sample_imputed2$SVTYPE, ]
    
    sample_truth3 = detected_SV[(GT_detected_SV[,i]%in% c("1/1","1|1")),]
    sample_imputed3 = imputed_SV[GT_imputed_SV[,i]=="1|1",]
    sample_TP3 <- sample_truth3[(sample_truth3$POS %in% as.vector(sapply(sample_imputed3$POS, function(x) (x-seq(-200,200))))) &
                                  sample_truth3$CHROM %in% sample_imputed3$CHROM &
                                  sample_truth3$SVTYPE %in% sample_imputed3$SVTYPE, ]
    
    TP <- nrow(sample_TP1) + nrow(sample_TP2) + nrow(sample_TP3)
    n_total_detected <- nrow(sample_truth1) + nrow(sample_truth2) + nrow(sample_truth3)
    recall <- TP/n_total_detected
    
    sample_TP1 <- sample_imputed1[(sample_imputed1$POS %in% as.vector(sapply(sample_truth1$POS, function(x) (x-seq(-200,200))))) &
                                    sample_imputed1$CHROM %in% sample_truth1$CHROM &
                                    sample_imputed1$SVTYPE %in% sample_truth1$SVTYPE,]
    sample_TP2 <- sample_imputed2[(sample_imputed2$POS %in% as.vector(sapply(sample_truth2$POS, function(x) (x-seq(-200,200))))) &
                                    sample_imputed2$CHROM %in% sample_truth2$CHROM&
                                    sample_imputed2$SVTYPE %in% sample_truth2$SVTYPE,]
    sample_TP3 <- sample_imputed3[(sample_imputed3$POS %in% as.vector(sapply(sample_truth3$POS, function(x) (x-seq(-200,200))))) &
                                    sample_imputed3$CHROM %in% sample_truth3$CHROM&
                                    sample_imputed3$SVTYPE %in% sample_truth3$SVTYPE,]
    
    TP <- nrow(sample_TP1) + nrow(sample_TP2) + nrow(sample_TP3)
    n_total_imputed <- nrow(sample_imputed1) + nrow(sample_imputed2) + nrow(sample_imputed3)
    precision <- TP/n_total_imputed
    
    df_SV <- rbind(df_SV, data.frame(TP, 
                                     truth = n_total_detected,
                                     imputed = n_total_imputed,
                                     recall, precision, 
                                     sample = truth_ONTsamples[i,1],
                                     var_type = "SV"))
  }
  
  return(df_SV)
}

imputation_TP_concordant_GT <- function(df_results,ref_name, threshold, NGS_detected, ONT_detected, input_name, DIR_sample_name){
  df_imputed <- df_results[df_results$IMP!="." & df_results$DR2>=threshold,]
  if(input_name == ""){
    sample_name <- ""
    N = 101
  }
  if(input_name == "CH_N100"){
    sample_name <- "_CH"
    N = 101
  }
  if(input_name == "CH_N60"){
    sample_name <- "_N60_CH"
    N = 61
  }
  if(input_name == "SAS_N60"){
    sample_name <- "_N60_SAS"
    N = 61
  }
  if(input_name == "EUR_N60"){
    sample_name <- "_N60_EUR"
    N = 61
  }
  if(input_name == "AFR_N60"){
    sample_name <- "_N60_AFR"
    N = 61
  }
  if(input_name == "AMR_N60"){
    sample_name <- "_N60_AMR"
    N = 61
  }
  
  if(input_name == "5groups_N150"){
    sample_name <- "_N150_5groups"
    N = 151
  }
  
  truth_ONTsamples <- read.table(paste0(DIR_sample_name, "truth_samples",sample_name,".txt"), header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  
  GT_imputed <- data.frame(stringr::str_split_fixed(df_imputed$GT," ",N))
  GT_imputed <- GT_imputed[,c(2:N)]
  colnames(GT_imputed) <- truth_ONTsamples[,1]
  
  GT_detected <- data.frame(stringr::str_split_fixed(NGS_detected$FORMAT_GT," ",N))
  GT_detected <- GT_detected[,c(2:N)]
  colnames(GT_detected) <- truth_ONTsamples[,1]
  
  df_imputed <- assign_newID_imputed(df_imputed, ref_name)
  df_imputed$MAF <- ifelse(df_imputed$AF>0.5, 1-df_imputed$AF, df_imputed$AF)
  
  imputed_ALL <- df_imputed
  imputed_SNV <- imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & (nchar(imputed_ALL$REF)==1) & (nchar(imputed_ALL$ALT)==1),]
  imputed_indel <- imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & !(nchar(imputed_ALL$REF)==1 & nchar(imputed_ALL$ALT)==1),]
  imputed_SV <- imputed_ALL[grepl("SV", imputed_ALL$ID),]
  
  GT_imputed_ALL <- GT_imputed
  GT_imputed_SNV <- GT_imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & (nchar(imputed_ALL$REF)==1) & (nchar(imputed_ALL$ALT)==1),]
  GT_imputed_indel <- GT_imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & !(nchar(imputed_ALL$REF)==1 & nchar(imputed_ALL$ALT)==1),]
  GT_imputed_SV <- GT_imputed_ALL[grepl("SV", imputed_ALL$ID),]
  
  detected_SNV <- NGS_detected[(nchar(NGS_detected$REF)==1) & (nchar(NGS_detected$ALT)==1), c(1:9)]
  detected_indel <- NGS_detected[!(nchar(NGS_detected$REF)==1 & nchar(NGS_detected$ALT)==1), c(1:9)]
  GT_detected_SNV <- GT_detected[(nchar(NGS_detected$REF)==1) & (nchar(NGS_detected$ALT)==1),]
  GT_detected_indel <- GT_detected[!(nchar(NGS_detected$REF)==1 & nchar(NGS_detected$ALT)==1),]
  
  detected_SV <- ONT_detected
  GT_detected_SV <- data.frame(stringr::str_split_fixed(ONT_detected$GT," ",N))
  GT_detected_SV <- GT_detected_SV[,c(2:N)]
  colnames(GT_detected_SV) <- truth_ONTsamples[,1]
  
  df_SNV <- imputation_TP_SNV_concordant_GT(detected_SNV, GT_detected_SNV, imputed_SNV, GT_imputed_SNV, truth_ONTsamples)
  df_indel <- imputation_TP_indel_concordant_GT(detected_indel, GT_detected_indel, imputed_indel, GT_imputed_indel, truth_ONTsamples)
  df_SV <- imputation_TP_SV_concordant_GT(detected_SV, GT_detected_SV, imputed_SV, GT_imputed_SV, truth_ONTsamples)
  
  df <- rbind(df_SNV, df_indel, df_SV)
  df$F1_score <- (2*df$recall*df$precision)/(df$recall+df$precision)
  
  return(df)
}

imputation_eval <- function(df_results, ref_name, threshold, NGS_detected, ONT_detected, input_name, GT_concordant, DIR_sample_name){
  if(GT_concordant){
    df <- imputation_TP_concordant_GT(df_results, ref_name, threshold, NGS_detected, ONT_detected, input_name, DIR_sample_name)
  }else{
    df <- imputation_TP(df_results, ref_name, threshold, NGS_detected, ONT_detected, input_name)
  }
  
  df2 <- data.frame(value = c(mean(df[df$var_type=="SNV",]$TP),
                              mean(df[df$var_type=="indel",]$TP),
                              mean(df[df$var_type=="SV",]$TP), 
                              mean(df[df$var_type=="SNV",]$truth),
                              mean(df[df$var_type=="indel",]$truth),
                              mean(df[df$var_type=="SV",]$truth), 
                              mean(df[df$var_type=="SNV",]$imputed),
                              mean(df[df$var_type=="indel",]$imputed),
                              mean(df[df$var_type=="SV",]$imputed), 
                              mean(df[df$var_type=="SNV",]$recall),
                              mean(df[df$var_type=="indel",]$recall),
                              mean(df[df$var_type=="SV",]$recall), 
                              mean(df[df$var_type=="SNV"& (!is.na(df$precision)),]$precision),
                              mean(df[df$var_type=="indel"& (!is.na(df$precision)),]$precision),
                              mean(df[df$var_type=="SV"& (!is.na(df$precision)),]$precision)),
                    
                    type = rep(c("TP","truth","imputed","recall", "precision"), each=3),
                    var_type = rep(c("SNV","indel","SV"), 5))
  
  df_eval <- df2
  return(df_eval)
}


imputation_eval_ALL <- function(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name, high_T, low_T){
  ref$MAF <- ifelse(ref$AF>0.5, 1-ref$AF, ref$AF)
  if(!("MAF" %in% colnames(SNV_truth))){
    SNV_truth$MAF <- ifelse(SNV_truth$AF>0.5, 1-SNV_truth$AF, SNV_truth$AF)
  }
  if(!("MAF" %in% colnames(SV_truth))){
    SV_truth$MAF <- ifelse(SV_truth$AF>0.5, 1-SV_truth$AF, SV_truth$AF)
  }
  
  tmp1 <- imputation_eval(df_results = ref[ref$MAF < high_T & ref$MAF >= low_T,], 
                          ref_name = ref_name, 
                          threshold = threshold, 
                          NGS_detected = SNV_truth, 
                          ONT_detected = SV_truth, 
                          input_name, GT_concordant, DIR_sample_name)
  tmp1[tmp1$type=="TP",]$type <- "TP_imputed"
  tmp2 <- imputation_eval(df_results = ref, 
                          ref_name = ref_name, 
                          threshold = threshold, 
                          NGS_detected = SNV_truth[SNV_truth$MAF>=low_T & SNV_truth$MAF<high_T,], 
                          ONT_detected = SV_truth[SV_truth$MAF>=low_T & SV_truth$MAF<high_T,],
                          input_name, GT_concordant, DIR_sample_name)
  tmp2[tmp2$type=="TP",]$type <- "TP_truth"
  
  df2 <- rbind(tmp1[tmp1$type %in% c("TP_imputed","imputed","precision"),],
               tmp2[tmp2$type %in% c("TP_truth","truth","recall"),] )
  
  SNV_F <- 2*df2[df2$var_type == "SNV" & df2$type == "recall",]$value*df2[df2$var_type == "SNV"& df2$type == "precision",]$value/
    (df2[df2$var_type == "SNV"& df2$type == "recall",]$value+df2[df2$var_type == "SNV"& df2$type == "precision",]$value)
  indel_F <-  2*df2[df2$var_type == "indel" & df2$type == "recall",]$value*df2[df2$var_type == "indel"& df2$type == "precision",]$value/
    (df2[df2$var_type == "indel"& df2$type == "recall",]$value+df2[df2$var_type == "indel"& df2$type == "precision",]$value)
  SV_F <- 2*df2[df2$var_type == "SV" & df2$type == "recall",]$value*df2[df2$var_type == "SV"& df2$type == "precision",]$value/
    (df2[df2$var_type == "SV"& df2$type == "recall",]$value+df2[df2$var_type == "SV"& df2$type == "precision",]$value)
  df_eval <- rbind(df2, 
                   data.frame(value = c(SNV_F, indel_F, SV_F),  
                              type = rep("F1_score", each=3),
                              var_type = c("SNV","indel","SV")))
  return(df_eval)
}


imputation_eval_byMAF <- function(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name){
  df_ALL <- imputation_eval_ALL(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name, high_T=0.51, low_T=0)
  df_MAF0.51 <- imputation_eval_ALL(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name, high_T=0.51, low_T=0.05)
  df_MAF0.05 <- imputation_eval_ALL(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name, high_T=0.05, low_T=0.01)
  df_MAF0.01 <- imputation_eval_ALL(ref, SNV_truth, SV_truth, ref_name, threshold, input_name, GT_concordant, DIR_sample_name, high_T=0.01, low_T=0)
  df <- cbind(rbind(cbind(df_ALL, ref=ref_name, freq = "ALL"),
                    cbind(df_MAF0.51, ref=ref_name, freq = "MAF>=0.05"),
                    cbind(df_MAF0.05, ref=ref_name, freq = "0.01<=MAF<0.05"),
                    cbind(df_MAF0.01, ref=ref_name, freq = "MAF<0.01")),
              DR2= paste0("DR2>=",threshold))
  return(df)
}

assign_newID_imputed <- function(df_imputed, ref_name){
  imputed_SV <- df_imputed[grepl("SV", df_imputed$ID),]
  
  if(ref_name == "ZYref"){
    imputed_SV$SVTYPE <- data.frame(stringr::str_split_fixed(imputed_SV$ID,"_",4))[,4]
  }
  if(ref_name == "1KGP"){
    imputed_SV$SVTYPE <- imputed_SV$ALT
    ref_ALT <- c("<DEL:ME", "<DEL>", "<INS:ME", "<INS>", "<INV>", "<DUP>")
    #std_svtype <- c("DEL_ME", "DEL", "INS_ME", "INS", "INV", "DUP")
    std_svtype <- c("DEL", "DEL", "INS", "INS", "INV", "DUP")
    for(i in c(1: length(ref_ALT))){
      if(sum(grepl(ref_ALT[i], imputed_SV$ALT)!=0)){
        imputed_SV[grepl(ref_ALT[i], imputed_SV$ALT),]$SVTYPE <- std_svtype[i]
      }
    }
  }
  if(ref_name =="ChinaMAP"){
    df_imputed$SVTYPE <- "."
    df_imputed$new_ID <- paste0(df_imputed$CHROM,"_", df_imputed$POS,"_",df_imputed$REF,"_", df_imputed$ALT)
  }else{
    df_imputed$SVTYPE <- "."
    df_imputed[grepl("SV", df_imputed$ID),]$SVTYPE <- imputed_SV$SVTYPE
    df_imputed$new_ID <- paste0(df_imputed$CHROM,"_", df_imputed$POS,"_",df_imputed$REF,"_", df_imputed$ALT)
    df_imputed[grepl("SV", df_imputed$ID),]$new_ID <- paste0(df_imputed[grepl("SV", df_imputed$ID),]$CHROM,
                                                             "_",df_imputed[grepl("SV", df_imputed$ID),]$POS,
                                                             "_SV_",df_imputed[grepl("SV", df_imputed$ID),]$SVTYPE)
  }
  return(df_imputed)
}
