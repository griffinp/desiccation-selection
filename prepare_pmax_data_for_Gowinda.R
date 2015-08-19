# The aim is to take the pmax data from Murali et al. 2014
# and convert it to a list of 'functional categories'
# in the form that Gowinda wants for GO analysis.

bodypart <- read.table("PERCENTMAX_FLYATLAS.txt", header=TRUE,
                       sep="\t", stringsAsFactors=FALSE)

dev_stage <- read.table("PERCENTMAX_MODENCODEFPKM.txt", header=TRUE,
                        sep="\t", stringsAsFactors=FALSE)

convert_pmax_table_to_gowinda_input_format <- function(pmax_table, filter_level){
  output_mat <- matrix(ncol=3)
  for(i in 2:ncol(pmax_table)){
    tempcolname <- colnames(pmax_table)[i]
    tempcol <- pmax_table[,i]
    tempfilter <- pmax_table[tempcol>=filter_level,1]
    gene_field <- paste(tempfilter, collapse=" ")
    output_mat<- rbind(output_mat, c(tempcolname,tempcolname,
                          gene_field))
  }
  return(as.data.frame(output_mat[2:nrow(output_mat),]))
}

bodypart_75pc <- convert_pmax_table_to_gowinda_input_format(bodypart, 75.0)
write.table(bodypart_75pc, "Pmax_bodypart_75pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

dev_stage_75pc <- convert_pmax_table_to_gowinda_input_format(dev_stage, 75.0)
write.table(dev_stage_75pc, "Pmax_devstage_75pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

bodypart_90pc <- convert_pmax_table_to_gowinda_input_format(bodypart, 90.0)
write.table(bodypart_90pc, "Pmax_bodypart_90pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

dev_stage_90pc <- convert_pmax_table_to_gowinda_input_format(dev_stage, 90.0)
write.table(dev_stage_90pc, "Pmax_devstage_90pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

bodypart_50pc <- convert_pmax_table_to_gowinda_input_format(bodypart, 50.0)
write.table(bodypart_50pc, "Pmax_bodypart_50pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

dev_stage_50pc <- convert_pmax_table_to_gowinda_input_format(dev_stage, 50.0)
write.table(dev_stage_50pc, "Pmax_devstage_50pc_filter.txt", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)





