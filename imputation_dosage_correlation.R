#!/usr/bin/env Rscript
#

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


imputation_dosage_SVs <- function(df_results,ref_name, threshold, DS_ONT_detected, input_name, DIR_sample_name){
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
  
  ##############Total Imputated dosage
  DS_imputed <- data.frame(stringr::str_split_fixed(df_imputed$DS," ",N))
  DS_imputed <- DS_imputed[,c(2:N)]
  colnames(DS_imputed) <- truth_ONTsamples[,1]
  
  ##############detected dosage SV
  colnames(DS_ONT_detected) <- c("CHROM","POS", truth_ONTsamples[,1],"ID")
  
  ####################Assign new ID
  df_imputed <- assign_newID_imputed(df_imputed, ref_name)
  df_imputed$MAF <- ifelse(df_imputed$AF>0.5, 1-df_imputed$AF, df_imputed$AF)
  DS_imputed$ID <- df_imputed$new_ID
  
  ####Imputed dosage
  DS_imputed_SV <- DS_imputed[grepl("SV", DS_imputed$ID),]
  DS_imputed_SV$CHROM <-  data.frame(stringr::str_split_fixed(DS_imputed_SV$ID,"_",4))[,1]
  DS_imputed_SV$POS <-  data.frame(stringr::str_split_fixed(DS_imputed_SV$ID,"_",4))[,2]
  DS_imputed_SV$SVTYPE <-  data.frame(stringr::str_split_fixed(DS_imputed_SV$ID,"_",4))[,4]
  
  #####Detected dosage
  DS_detected_SV <- DS_ONT_detected
  DS_detected_SV$SVTYPE <- data.frame(stringr::str_split_fixed(DS_detected_SV$ID,"_",4))[,4]
  
  ######calcuate for SVs
  dosage_rho_SV <- c()
  matched_imputed_SV <- SVmatch_detected_imputed(DS_detected_SV, DS_imputed_SV)
  index <- match(DS_detected_SV$ID, matched_imputed_SV$new_ID)
  
  for(i in c(1:nrow(truth_ONTsamples))){
    x = DS_detected_SV[!is.na(index),][,i+2]
    y = as.numeric(matched_imputed_SV[index[!is.na(index)],][,i])
    result <- cor(x,y, method = "pearson")
    dosage_rho_SV <- rbind(dosage_rho_SV, 
                           data.frame(r_squared=(result)^2, 
                                      sampleID = truth_ONTsamples[i,1],
                                      vaiant_assessed = length(x)))
  }
  dosage_rho <- cbind(dosage_rho_SV, variant_type = "SV")
  return(dosage_rho)
}


SVmatch_detected_imputed <- function(DS_detected_SV, DS_imputed_SV){
  matched_imputed_SV <- c()
  for(i in c(1:nrow(DS_detected_SV))){
    tmp <- DS_imputed_SV[(DS_imputed_SV$POS %in% as.vector(sapply(as.numeric(DS_detected_SV$POS[i]), function(x) (x-seq(-200,200))))) &
                           DS_imputed_SV$CHROM == DS_detected_SV$CHROM[i] &
                           DS_imputed_SV$SVTYPE == DS_detected_SV$SVTYPE[i],]
    if(nrow(tmp)!=0){
      tmp$new_ID <- DS_detected_SV$ID[i]
      matched_imputed_SV <- rbind(matched_imputed_SV, tmp)
    }
  }
  return(matched_imputed_SV)
}

imputation_dosage_smallvariants_byAF <- function(df_results, ref_name, threshold, DS_NGS_detected, input_name, DIR_sample_name){
  ALL <- df_results
  high_T <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  low_T <- c(0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  df_AF <- c()
  for(i in c(1:9)){
    df <- imputation_dosage_smallvariants(df_results = ALL[ALL$AF>=low_T[i] & ALL$AF<high_T[i],], 
                                          ref_name = ref_name, 
                                          threshold = threshold, 
                                          DS_NGS_detected,
                                          input_name, 
                                          DIR_sample_name)
    df_AF <- rbind(df_AF, cbind(df,low_T = low_T[i], high_T = high_T[i]))
  }
  return(df_AF)
}
imputation_dosage_smallvariants_byMAF <- function(df_results, ref_name, threshold, DS_NGS_detected,input_name, DIR_sample_name){
  ALL <- df_results
  ALL$MAF <-  ifelse(ALL$AF>0.5, 1-ALL$AF, ALL$AF)
  high_T <- c(0.01, 0.05, 0.51)
  low_T <- c(0, 0.01,0.05)
  df_AF <- c()
  for(i in c(1:3)){
    df <- imputation_dosage_smallvariants(df_results = ALL[ALL$MAF>=low_T[i] & ALL$MAF<high_T[i],], 
                                          ref_name = ref_name, 
                                          threshold = threshold, 
                                          DS_NGS_detected,
                                          input_name, 
                                          DIR_sample_name)
    df_AF <- rbind(df_AF, cbind(df,low_T = low_T[i], high_T = high_T[i]))
  }
  return(df_AF)
}

imputation_dosage_SVs_byAF <- function(df_results, ref_name, threshold, DS_ONT_detected, input_name, DIR_sample_name){
  ALL <- df_results
  high_T <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  low_T <- c(0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  df_AF <- c()
  for(i in c(1:9)){
    df <- imputation_dosage_SVs(df_results = ALL[ALL$AF>=low_T[i] & ALL$AF<high_T[i],], 
                                ref_name = ref_name, 
                                threshold = threshold, 
                                DS_ONT_detected, 
                                input_name, 
                                DIR_sample_name)
    df_AF <- rbind(df_AF, cbind(df,low_T = low_T[i], high_T = high_T[i]))
  }
  return(df_AF)
}

imputation_dosage_SVs_byMAF <- function(df_results, ref_name, threshold, DS_ONT_detected, input_name, DIR_sample_name){
  ALL <- df_results
  ALL$MAF <-  ifelse(ALL$AF>0.5, 1-ALL$AF, ALL$AF)
  high_T <- c(0.01, 0.05, 0.51)
  low_T <- c(0, 0.01,0.05)
  df_AF <- c()
  for(i in c(1:3)){
    df <- imputation_dosage_SVs(df_results = ALL[ALL$MAF>=low_T[i] & ALL$MAF<high_T[i],], 
                                ref_name = ref_name, 
                                threshold = threshold, 
                                DS_ONT_detected, 
                                input_name, 
                                DIR_sample_name)
    df_AF <- rbind(df_AF, cbind(df,low_T = low_T[i], high_T = high_T[i]))
  }
  return(df_AF)
}


imputation_dosage_smallvariants <- function(df_results,ref_name, threshold, DS_NGS_detected, input_name, DIR_sample_name){
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
  
  ##############Total Imputated dosage
  DS_imputed <- data.frame(stringr::str_split_fixed(df_imputed$DS," ",N))
  DS_imputed <- DS_imputed[,c(2:N)]
  colnames(DS_imputed) <- truth_ONTsamples[,1]
  
  ##############Total detected dosage small variants
  colnames(DS_NGS_detected) <- c("CHROM","POS","REF","ALT", truth_ONTsamples[,1])
  DS_NGS_detected$ID <- paste0(DS_NGS_detected$CHROM,"_", DS_NGS_detected$POS,"_",DS_NGS_detected$REF,"_", DS_NGS_detected$ALT)
  
  ####################Assign new ID
  df_imputed <- assign_newID_imputed(df_imputed, ref_name)
  df_imputed$MAF <- ifelse(df_imputed$AF>0.5, 1-df_imputed$AF, df_imputed$AF)
  DS_imputed$ID <- df_imputed$new_ID
  DS_imputed$REF <- df_imputed$REF
  DS_imputed$ALT <- df_imputed$ALT
  
  ####Imputed dosage
  DS_imputed_SNV <- DS_imputed[(!grepl("SV", DS_imputed$ID)) & (nchar(DS_imputed$REF)==1) & (nchar(DS_imputed$ALT)==1),]
  DS_imputed_indel <- DS_imputed[(!grepl("SV", DS_imputed$ID)) & !(nchar(DS_imputed$REF)==1 & nchar(DS_imputed$ALT)==1),]
  
  #####Detected dosage
  DS_detected_SNV <- DS_NGS_detected[(nchar(DS_NGS_detected$REF)==1) & (nchar(DS_NGS_detected$ALT)==1),]
  DS_detected_indel <- DS_NGS_detected[!(nchar(DS_NGS_detected$REF)==1 & nchar(DS_NGS_detected$ALT)==1),]
  
  ######calculate for SNVs
  dosage_rho_SNV <- c()
  index <- match(DS_detected_SNV$ID, DS_imputed_SNV$ID)
  matched_detected <- DS_detected_SNV[!is.na(index),]
  matched_imputed <- DS_imputed_SNV[index[!is.na(index)],]
  
  for(i in c(1:nrow(truth_ONTsamples))){
    x = as.numeric(matched_detected[, (i+4)])
    y = as.numeric(matched_imputed[, i])
    result <- cor.test(x,y, method = "pearson")
    dosage_rho_SNV <- rbind(dosage_rho_SNV, 
                            data.frame(r_squared=(result$estimate)^2, 
                                       sampleID = truth_ONTsamples[i,1],
                                       vaiant_assessed = length(x)))
  }
  
  ######calculate for indels
  if(ref_name == "ChinaMAP"){
    dosage_rho <- cbind(dosage_rho_SNV, variant_type = "SNV")
  }else{
    dosage_rho_indel <- c()
    index <- match(DS_detected_indel$ID, DS_imputed_indel$ID)
    matched_detected <- DS_detected_indel[!is.na(index),]
    matched_imputed <- DS_imputed_indel[index[!is.na(index)],]
    
    for(i in c(1:nrow(truth_ONTsamples))){
      x = as.numeric(matched_detected[, (i+4)])
      y = as.numeric(matched_imputed[, i])
      result <- cor.test(x,y, method = "pearson")
      dosage_rho_indel <- rbind(dosage_rho_indel, 
                                data.frame(r_squared=(result$estimate)^2, 
                                           sampleID = truth_ONTsamples[i,1],
                                           vaiant_assessed = length(x)))
    }
    
    dosage_rho <- rbind(cbind(dosage_rho_SNV, variant_type = "SNV"),
                        cbind(dosage_rho_indel, variant_type = "indel"))
  }
  
  return(dosage_rho)
}
