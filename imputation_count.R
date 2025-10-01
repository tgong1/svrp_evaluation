#!/usr/bin/env Rscript
#
generage_imputation_count <- function(imputed_ALL, var, ZY_ref){
  imputed_SNV <- imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & (nchar(imputed_ALL$REF)==1) & (nchar(imputed_ALL$ALT)==1),]
  imputed_indel <- imputed_ALL[(!grepl("SV", imputed_ALL$ID)) & !(nchar(imputed_ALL$REF)==1 & nchar(imputed_ALL$ALT)==1),]
  imputed_SV <- imputed_ALL[grepl("SV", imputed_ALL$ID),]
  
  if(nrow(imputed_SV)==0){
    df_svtype <- c()
  }else{
    if(ZY_ref == TRUE){
      tmp <- data.frame(stringr::str_split_fixed(imputed_SV$ID,"_",4))[,4]
      imputed_SV$SVTYPE <- sub(x=tmp,"_.*","")
    }else{
      imputed_SV$SVTYPE <- imputed_SV$ALT
      ref_ALT <- c("<DEL:ME", "<DEL>", "<INS:ME", "<INS>", "<INV>", "<DUP>")
      std_svtype <- c("DEL_ME", "DEL", "INS_ME", "INS", "INV", "DUP")
      for(i in c(1: length(ref_ALT))){
        if(sum(grepl(ref_ALT[i], imputed_SV$ALT)!=0)){
          imputed_SV[grepl(ref_ALT[i], imputed_SV$ALT),]$SVTYPE <- std_svtype[i]
        }
      }
    }
    
    df_svtype <- data.frame(table(imputed_SV$SVTYPE))
    colnames(df_svtype) <- c("var_type","count")
    df_svtype$var_type <- as.character(df_svtype$var_type)
  }
  df_imputed <- data.frame(rbind(c("total", nrow(imputed_ALL)),
                                 c("SNV", nrow(imputed_SNV)),
                                 c("indel", nrow(imputed_indel)),
                                 c("SV", nrow(imputed_SV)), 
                                 df_svtype))
  df_imputed$var <- var
  colnames(df_imputed) <- c("var_type","count","var")
  return(df_imputed)
}

imputation_count_target <- function(df_results, ZY_ref, T_low, T_high, AF_T){
  colnames(df_results) <- c("CHROM","POS","ID","REF","ALT","DR2","AF","IMP")
  
  if(AF_T){
    df_results <- df_results[df_results$AF>0,]
  }
  df_results$AF <- as.numeric(df_results$AF)
  df_results$MAF <- ifelse(df_results$AF>0.5, 1-df_results$AF, df_results$AF)
  
  df_results <- df_results[df_results$IMP!=".",]
  df_results$DR2 <- as.numeric(df_results$DR2)
  
  imputed_ALL <- df_results[df_results$DR2>=T_low & df_results$DR2< T_high,]
  df_imputed <- rbind(generage_imputation_count(imputed_ALL, var = "all", ZY_ref),
                      generage_imputation_count(imputed_ALL[imputed_ALL$MAF >= 0.05,], var = "MAF>=0.05", ZY_ref),
                      generage_imputation_count(imputed_ALL[imputed_ALL$MAF >= 0.01 & imputed_ALL$MAF < 0.05,], var = "0.01<=MAF<0.05", ZY_ref),
                      generage_imputation_count(imputed_ALL[imputed_ALL$MAF < 0.01,], var = "MAF<0.01", ZY_ref))
  df_imputed$count <- as.numeric(df_imputed$count)
  df_imputed$DR2 <- paste0(T_low,"<=DR2<", T_high)
  return(df_imputed)
}
