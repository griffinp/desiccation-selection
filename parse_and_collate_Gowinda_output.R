#################################################
# Reading in Gowinda GO annotation output files #
#################################################

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment/Gowinda_files")

select_sig_category <- function(Gowinda_output, sig_threshold, uncorrected_or_FDR){
  if(uncorrected_or_FDR=='uncorrected'){
    sig_category <- subset(Gowinda_output, Gowinda_output[,4]<sig_threshold)
    output_table <- sig_category[,c(9,4,5)]
    colnames(output_table) <- c("GO_category", "Uncorrected_pval", "FDR_pval")
    return(output_table)
  }
  if(uncorrected_or_FDR=='FDR'){
    sig_category <- subset(Gowinda_output, Gowinda_output[,5]<sig_threshold)
    output_table <- sig_category[,c(9,4,5)]
    colnames(output_table) <- c("GO_category", "Uncorrected_pval", "FDR_pval")
    return(output_table)
  }
}

collate_results <- function(input_object_list){
  for(i in 1:length(input_object_list)){
    temp_table <- input_object_list[[i]]
    temp_ordered <- temp_table[order(temp_table[,1]),]
    if(i==1){
      output_fdr <- data.frame(temp_ordered[,5])
      output_uncorrected <- data.frame(temp_ordered[,4])
      rownames(output_fdr) <- temp_ordered[,1]
      rownames(output_uncorrected) <- temp_ordered[,1]
    }
    else{
      output_fdr <- cbind(output_fdr, temp_ordered[,5])
      output_uncorrected <- cbind(output_uncorrected, temp_ordered[,4])
    }
  }
  output_table <- cbind(output_fdr, output_uncorrected)
  colnames(output_table) <- paste(rep(c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5"), times=2), 
                                  rep(c("FDR_pval", "Uncorrected_pval"), each=10), sep="_")
  return(output_table)
}

possible_devstage_levels <- c('ADULT_FEMALE_ECLOSION_PLUS_1D', 'ADULT_FEMALE_ECLOSION_PLUS_30D', 
                              'ADULT_FEMALE_ECLOSION_PLUS_5D', 'ADULT_MALE_ECLOSION_PLUS_1D', 
                              'ADULT_MALE_ECLOSION_PLUS_30D', 'ADULT_MALE_ECLOSION_PLUS_5D', 
                              'EMBRYOS_0_2_HR', 'EMBRYOS_10_12_HR', 'EMBRYOS_12_14_HR', 'EMBRYOS_14_16_HR', 
                              'EMBRYOS_16_18_HR', 'EMBRYOS_18_20_HR', 'EMBRYOS_2_4_HR', 'EMBRYOS_20_22_HR', 
                              'EMBRYOS_22_24_HR', 'EMBRYOS_4_6_HR', 'EMBRYOS_6_8_HR', 'EMBRYOS_8_10_HR', 
                              'L1_LARVAE', 'L2_LARVAE', 'L3_LARVAE_12_HR_POST_MOLT', 
                              'L3_LARVAE_CLEAR_GUT_PS_7_9', 'L3_LARVAE_DARK_BLUE_GUT_PS_1_2', 
                              'L3LARVAE_LIGHT_BLUE_GUT_PS_3_6', 'PUPAE_WPP_PLUS_2D', 'PUPAE_WPP_PLUS_3D', 
                              'PUPAE_WPP_PLUS_4D', 'WHITE_PRE_PUPAE', 'WPP_PLUS_12_HR', 'WPP_PLUS_24_HR')

possible_bodypart_levels <- c('ADULT_CARCASS', 'ADULT_FATBODY', 'ADULT_SALIVARY_GLAND', 'BRAIN', 
                              'CROP', 'EYE', 'HEART', 'HINDGUT', 'LARVAL_CARCASS', 'LARVAL_CNS', 
                              'LARVAL_FATBODY', 'LARVAL_HINDGUT', 'LARVAL_MIDGUT', 
                              'LARVAL_SALIVARY_GLAND', 'LARVAL_TRACHEA', 'LARVAL_TUBULE', 
                              'MALE_ACCESSORY_GLAND', 'MIDGUT', 'OVARY', 'SPERMATHECA_VIRGIN', 
                              'TESTIS', 'THORACICOABDOMINAL_GANGLION', 'TUBULE')



for(i in Sample_code){
  temp_go_name <- paste(i, "GO", sep="_")
  temp_bp_50_name <- paste(i, "bodypart_50pc", sep="_")
  temp_bp_75_name <- paste(i, "bodypart_75pc", sep="_")
  temp_bp_90_name <- paste(i, "bodypart_90pc", sep="_")
  temp_ds_50_name <- paste(i, "devstage_50pc", sep="_")
  temp_ds_75_name <- paste(i, "devstage_75pc", sep="_")
  temp_ds_90_name <- paste(i, "devstage_90pc", sep="_")
  if(i%in%c("C1", "C2", "C3", "C4", "C5")){
    temp_go_filename <- paste(i, "1000bp_up_down_output.txt", sep="_")
    temp_bp_50_filename <- paste(i, "1000bp_bodypart_50pc_output.txt", sep="_")
    temp_bp_75_filename <- paste(i, "1000bp_bodypart_75pc_output.txt", sep="_")
    temp_bp_90_filename <- paste(i, "1000bp_bodypart_90pc_output.txt", sep="_")
    temp_ds_50_filename <- paste(i, "1000bp_devstage_50pc_output.txt", sep="_")
    temp_ds_75_filename <- paste(i, "1000bp_devstage_75pc_output.txt", sep="_")
    temp_ds_90_filename <- paste(i, "1000bp_devstage_90pc_output.txt", sep="_")
  } else {
    temp_go_filename <- paste("noC", i, "1000bp_up_down_output.txt", sep="_")
    temp_bp_50_filename <- paste("noC", i, "1000bp_bodypart_50pc_output.txt", sep="_")
    temp_bp_75_filename <- paste("noC", i, "1000bp_bodypart_75pc_output.txt", sep="_")
    temp_bp_90_filename <- paste("noC", i, "1000bp_bodypart_90pc_output.txt", sep="_")
    temp_ds_50_filename <- paste("noC", i, "1000bp_devstage_50pc_output.txt", sep="_")
    temp_ds_75_filename <- paste("noC", i, "1000bp_devstage_75pc_output.txt", sep="_")
    temp_ds_90_filename <- paste("noC", i, "1000bp_devstage_90pc_output.txt", sep="_")
  }
  assign(temp_go_name, read.table(temp_go_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_bp_50_name, read.table(temp_bp_50_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_bp_75_name, read.table(temp_bp_75_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_bp_90_name, read.table(temp_bp_90_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_ds_50_name, read.table(temp_ds_50_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_ds_75_name, read.table(temp_ds_75_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  assign(temp_ds_90_name, read.table(temp_ds_90_filename, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE))
  
  temp_sig_GO <- select_sig_category(get(temp_go_name), 0.05, 'uncorrected')
  write.table(temp_sig_GO, file=paste(i, "sig_GO_category.txt", sep="_"), quote=FALSE, 
              sep="\t", row.names=FALSE)
  
  for(j in c(temp_bp_50_name, temp_bp_75_name, temp_bp_90_name)){
    j_df <- get(j)
    missing_levels <- possible_bodypart_levels[possible_bodypart_levels%in%j_df[,1]==FALSE]
    if(length(missing_levels)>0){
      missing_df <- as.data.frame(matrix(NA, ncol=10, nrow=length(missing_levels)))
      missing_df[,1] <- missing_levels
      j_df <- rbind(j_df, missing_df)
    }
    assign(j, j_df)
  }
  for(k in c(temp_ds_50_name, temp_ds_75_name, temp_ds_90_name)){
    k_df <- get(k)
    missing_levels <- possible_devstage_levels[possible_devstage_levels%in%k_df[,1]==FALSE]
    if(length(missing_levels)>0){
      missing_df <- as.data.frame(matrix(NA, ncol=10, nrow=length(missing_levels)))
      missing_df[,1] <- missing_levels
      k_df <- rbind(k_df, missing_df)
    }
    assign(k, k_df)
  }
}

#Once above code block has been run, all objects should be imported

bodypart_90pc_output <- collate_results(list(C1_bodypart_90pc, C2_bodypart_90pc,
                                             C3_bodypart_90pc, C4_bodypart_90pc,
                                             C5_bodypart_90pc, D1_bodypart_90pc, 
                                             D2_bodypart_90pc, D3_bodypart_90pc, 
                                             D4_bodypart_90pc, D5_bodypart_90pc))
write.table(bodypart_90pc_output, file="Bodypart_90pc_noC_D_collated.txt", quote=FALSE)

devstage_90pc_output <- collate_results(list(C1_devstage_90pc, C2_devstage_90pc,
                                             C3_devstage_90pc, C4_devstage_90pc,
                                             C5_devstage_90pc, D1_devstage_90pc, 
                                             D2_devstage_90pc, D3_devstage_90pc, 
                                             D4_devstage_90pc, D5_devstage_90pc))
write.table(devstage_90pc_output, file="devstage_90pc_noC_D_collated.txt", quote=FALSE)










# C1_GO <- read.table("C1_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C2_GO <- read.table("C2_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C3_GO <- read.table("C3_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C4_GO <- read.table("C4_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C5_GO <- read.table("C5_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D1_GO <- read.table("noC_D1_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D2_GO <- read.table("noC_D2_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D3_GO <- read.table("noC_D3_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D4_GO <- read.table("noC_D4_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D5_GO <- read.table("noC_D5_1000bp_up_down_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# 
# C1_bodypart_50pc <- read.table("C1_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C2_bodypart_50pc <- read.table("C2_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C3_bodypart_50pc <- read.table("C3_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C4_bodypart_50pc <- read.table("C4_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C5_bodypart_50pc <- read.table("C5_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D1_bodypart_50pc <- read.table("noC_D1_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D2_bodypart_50pc <- read.table("noC_D2_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D3_bodypart_50pc <- read.table("noC_D3_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D4_bodypart_50pc <- read.table("noC_D4_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D5_bodypart_50pc <- read.table("noC_D5_1000bp_bodypart_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C1_devstage_50pc <- read.table("C1_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C2_devstage_50pc <- read.table("C2_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C3_devstage_50pc <- read.table("C3_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C4_devstage_50pc <- read.table("C4_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# C5_devstage_50pc <- read.table("C5_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D1_devstage_50pc <- read.table("noC_D1_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D2_devstage_50pc <- read.table("noC_D2_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D3_devstage_50pc <- read.table("noC_D3_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D4_devstage_50pc <- read.table("noC_D4_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
# D5_devstage_50pc <- read.table("noC_D5_1000bp_devstage_50pc_output.txt", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)

# C1_sig_GO <- select_sig_GO(C1_GO, 0.05, 'uncorrected')
# write.table(C1_sig_GO, file="C1_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# C2_sig_GO <- select_sig_GO(C2_GO, 0.05, 'uncorrected')
# write.table(C2_sig_GO, file="C2_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# C3_sig_GO <- select_sig_GO(C3_GO, 0.05, 'uncorrected')
# write.table(C3_sig_GO, file="C3_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# C4_sig_GO <- select_sig_GO(C4_GO, 0.05, 'uncorrected')
# write.table(C4_sig_GO, file="C4_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# C5_sig_GO <- select_sig_GO(C5_GO, 0.05, 'uncorrected')
# write.table(C5_sig_GO, file="C5_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# 
# D1_sig_GO <- select_sig_GO(D1_GO, 0.05, 'uncorrected')
# write.table(D1_sig_GO, file="D1_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# D2_sig_GO <- select_sig_GO(D2_GO, 0.05, 'uncorrected')
# write.table(D2_sig_GO, file="D2_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# D3_sig_GO <- select_sig_GO(D3_GO, 0.05, 'uncorrected')
# write.table(D3_sig_GO, file="D3_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# D4_sig_GO <- select_sig_GO(D4_GO, 0.05, 'uncorrected')
# write.table(D4_sig_GO, file="D4_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)
# D5_sig_GO <- select_sig_GO(D5_GO, 0.05, 'uncorrected')
# write.table(D5_sig_GO, file="D5_sig_GO_categories.txt", quote=FALSE, 
#             sep="\t", row.names=FALSE)





