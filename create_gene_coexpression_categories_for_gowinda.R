# This script aims to make a 'gene category' file for use in Gowinda
# using the modules of co-expressed genes reported in Huang et al. (2015) PNAS

setwd("~/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment")

female_modules <- read.csv("Huang_et_al_female_modules.csv", header=TRUE,
                           sep=",", stringsAsFactors=FALSE)

module_length <- tapply(X=female_modules$FlyBase.ID, INDEX=female_modules$Module, 
                        FUN=length)
# long_modules <- names(which(module_length>=10))
# 
# female_long_modules <- subset(female_modules, 
#                               as.character(female_modules$Module)%in%long_modules)

gowinda_format <- data.frame(category_name=character(), 
                             category_name_2=character(),
                             gene_members=character())
for(i in 1:length(module_length)){
  temp_module <- names(module_length)[i]
  temp_module_name <- paste("female_module_", temp_module, sep="")
  temp_FBgn_names <- female_modules[which(female_modules$Module==as.integer(temp_module)), 1]
  temp_FBgn_as_vector <- paste(temp_FBgn_names, collapse=" ")
  temp_data_frame <- data.frame(category_name=temp_module_name, 
                                category_name_2=temp_module_name,
                                gene_members=temp_FBgn_as_vector)
  gowinda_format <- rbind(gowinda_format, temp_data_frame)
}

write.table(gowinda_format, file="female_coexpression_modules.txt", col.names=FALSE,
            row.names=FALSE, sep="\t", quote=FALSE)



male_modules <- read.csv("Huang_et_al_male_modules.csv", header=TRUE,
                           sep=",", stringsAsFactors=FALSE)

male_module_length <- tapply(X=male_modules$FlyBase.ID, INDEX=male_modules$Module, 
                        FUN=length)

male_gowinda_format <- data.frame(category_name=character(), 
                             category_name_2=character(),
                             gene_members=character())
for(i in 1:length(male_module_length)){
  temp_module <- names(male_module_length)[i]
  temp_module_name <- paste("male_module_", temp_module, sep="")
  temp_FBgn_names <- male_modules[which(male_modules$Module==as.integer(temp_module)), 1]
  temp_FBgn_as_vector <- paste(temp_FBgn_names, collapse=" ")
  temp_data_frame <- data.frame(category_name=temp_module_name, 
                                category_name_2=temp_module_name,
                                gene_members=temp_FBgn_as_vector)
  male_gowinda_format <- rbind(male_gowinda_format, temp_data_frame)
}

both_sexes <- rbind(gowinda_format, male_gowinda_format)

write.table(male_gowinda_format, file="male_coexpression_modules.txt", col.names=FALSE,
            row.names=FALSE, sep="\t", quote=FALSE)

write.table(both_sexes, file="both_sexes_coexpression_modules.txt", col.names=FALSE,
            row.names=FALSE, sep="\t", quote=FALSE)

