# Extract the genes that overlap among D1, D4, D5 and test for functional enrichment

########
# MAIN #
########

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

# Get object containing simulated gene lists and overlap lists
# (this is the version saved *after* parsing 'exon' gene names)

exons_parsed <- readRDS(file="/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/Gene_lists_exons_parsed_including_overlap_from_resampled_SNPs.rds")

# Extract overlap genes called in each combination of D replicates

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

# Obtain full list of sig C genes

C_sig_genes <- read.csv('Putative_lab_adaptation_gene_list.txt', header=FALSE,
                        stringsAsFactors=FALSE)[,1]

# Now get D1, D4, D5 gene lists

D1 <- read.table(file="D1_noC_sig_gene_names_FBgn_format.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE)
D4 <- read.table(file="D4_noC_sig_gene_names_FBgn_format.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE)
D5 <- read.table(file="D5_noC_sig_gene_names_FBgn_format.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE)

D1_D4<-intersect(D1[,1], D4[,1])
D1_D5<-intersect(D1[,1], D5[,1])
D4_D5<-intersect(D4[,1], D5[,1])
D1_D4_D5<-intersect(intersect(D1[,1], D4[,1]), D5[,1])

write.table(D1_D4_D5, file="D1_D4_D5_gene_list_overlap_FBgn.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Import Rahul's KEGG pathway databases (both 'empirical' and 'putative')

emp <- read.table("/Users/pgriffin/Documents/Drosophila Selection Experiment/rahul_gene_feature_enrichment/DMEL.KEGG.empirical.DB.tab", header=FALSE,
                  stringsAsFactors=FALSE)
put <- read.table("/Users/pgriffin/Documents/Drosophila Selection Experiment/rahul_gene_feature_enrichment/DMEL.KEGG.putative.DB.tab", header=FALSE,
                  stringsAsFactors=FALSE)

# Import body part category list

body <- read.table("/Users/pgriffin/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment/Gowinda_files/Pmax_bodypart_75pc_filter.txt", header=FALSE,
                   stringsAsFactors=FALSE, sep="\t")

ovary_genes <- unlist(strsplit(body[which(body[,1]=="OVARY"),3], split=" "))
D1_D4_ovary_genes <- D1_D4[which(D1_D4%in%ovary_genes)]

D1_D4_D5_ovary_genes <- D1_D4_D5[which(D1_D4_D5%in%ovary_genes)]

limonene_degradation <- put[put[,1]=="dme00903",]
lim_genes <- unlist(strsplit(limonene_degradation[,2], split=","))

galactose_metabolism <- put[put[,1]=="dme00052",]
gal_genes <- unlist(strsplit(galactose_metabolism[,2], split=","))

gal_genes[which(gal_genes%in%D1_D4_D5)]
