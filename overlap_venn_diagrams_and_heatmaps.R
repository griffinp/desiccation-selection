library(VennDiagram)

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

for(i in Sample_code){
  snp_temp_table_name <- paste(i, "sig_positions_table", sep="_")
  snp_temp_object_name <- paste(i, "sig_positions", sep="_")
  snp_temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
                          i, "_sig_SNPs.txt", sep="")
  assign(snp_temp_table_name, read.table(file=snp_temp_file_name, 
                                         header=FALSE, sep="\t",
                                         stringsAsFactors=FALSE))
  assign(snp_temp_object_name, paste(get(snp_temp_table_name)[,1], get(snp_temp_table_name)[,2], sep="_"))
  #setting no field separator so that chr and pos occupy one field
  
  #now processing gene name lists
  #using lists of exon-parsed genes, not FBgn format
  gene_temp_table_name <- paste(i, "sig_gene_names", sep="_")
  if(i %in% c("D1", "D2", "D3", "D4", "D5")){
    gene_temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
                                i, "_noC_sig_gene_list.txt", sep="")
  } else {
    gene_temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
                                 i, "_sig_gene_list.txt", sep="")
  }
  assign(gene_temp_table_name, read.table(file=gene_temp_file_name, 
                                          header=FALSE, quote="\"",
                                          stringsAsFactors=FALSE)[,1])  
  # Also importing gene lists for D replicates including C genes
  withC_table_name <- paste(i, "withC_gene_names", sep="_")
  if(i %in% c("D1", "D2", "D3", "D4", "D5")){
    withC_gene_temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
                                      i, "_sig_withC_gene_list.txt", sep="")
    assign(withC_table_name, read.table(file=withC_gene_temp_file_name, 
                                        header=FALSE, quote="\"",
                                        stringsAsFactors=FALSE)[,1])  
  }
}

############################
# SNP Venn Diagram objects #
############################

C_venn_SNPs <- draw.quintuple.venn(area1=length(C1_sig_positions), area2=length(C2_sig_positions), area3=length(C3_sig_positions),
                    area4=length(C4_sig_positions), area5=length(C5_sig_positions),
                 n12=length(intersect(C1_sig_positions, C2_sig_positions)), 
                 n13=length(intersect(C1_sig_positions, C3_sig_positions)),
                 n14=length(intersect(C1_sig_positions, C4_sig_positions)),
                 n15=length(intersect(C1_sig_positions, C5_sig_positions)),
                 n23=length(intersect(C2_sig_positions, C3_sig_positions)),
                 n24=length(intersect(C2_sig_positions, C4_sig_positions)),
                 n25=length(intersect(C2_sig_positions, C5_sig_positions)),
                 n34=length(intersect(C3_sig_positions, C4_sig_positions)),
                 n35=length(intersect(C3_sig_positions, C5_sig_positions)),
                 n45=length(intersect(C4_sig_positions, C5_sig_positions)),
                 n123=length(intersect(intersect(C1_sig_positions, C2_sig_positions), C3_sig_positions)),
                 n124=length(intersect(intersect(C1_sig_positions, C2_sig_positions), C4_sig_positions)),
                 n125=length(intersect(intersect(C1_sig_positions, C2_sig_positions), C5_sig_positions)),
                 n134=length(intersect(intersect(C1_sig_positions, C3_sig_positions), C4_sig_positions)),
                 n135=length(intersect(intersect(C1_sig_positions, C3_sig_positions), C5_sig_positions)),
                 n145=length(intersect(intersect(C1_sig_positions, C4_sig_positions), C5_sig_positions)),
                 n234=length(intersect(intersect(C2_sig_positions, C3_sig_positions), C4_sig_positions)),
                 n235=length(intersect(intersect(C2_sig_positions, C3_sig_positions), C5_sig_positions)),
                 n245=length(intersect(intersect(C2_sig_positions, C4_sig_positions), C5_sig_positions)),
                 n345=length(intersect(intersect(C3_sig_positions, C4_sig_positions), C5_sig_positions)),
                 n1234=length(intersect(intersect(C1_sig_positions, C2_sig_positions), 
                                        intersect(C3_sig_positions, C4_sig_positions))),
                 n1235=length(intersect(intersect(C1_sig_positions, C2_sig_positions),
                                        intersect(C3_sig_positions, C5_sig_positions))),
                 n1245=length(intersect(intersect(C1_sig_positions, C2_sig_positions),
                                        intersect(C4_sig_positions, C5_sig_positions))),
                 n1345=length(intersect(intersect(C1_sig_positions, C3_sig_positions),
                                        intersect(C4_sig_positions, C5_sig_positions))),
                 n2345=length(intersect(intersect(C2_sig_positions, C3_sig_positions),
                                        intersect(C4_sig_positions, C5_sig_positions))),
                 n12345=length(intersect(intersect(C1_sig_positions, C2_sig_positions), 
                                         intersect(intersect(C3_sig_positions, C4_sig_positions), C5_sig_positions))),
                 category=c("C1", "C2", "C3", "C4", "C5"))

 
D_venn_SNPs <- draw.quintuple.venn(area1=length(D1_sig_positions), area2=length(D2_sig_positions), area3=length(D3_sig_positions),
                                   area4=length(D4_sig_positions), area5=length(D5_sig_positions),
                                   n12=length(intersect(D1_sig_positions, D2_sig_positions)), 
                                   n13=length(intersect(D1_sig_positions, D3_sig_positions)),
                                   n14=length(intersect(D1_sig_positions, D4_sig_positions)),
                                   n15=length(intersect(D1_sig_positions, D5_sig_positions)),
                                   n23=length(intersect(D2_sig_positions, D3_sig_positions)),
                                   n24=length(intersect(D2_sig_positions, D4_sig_positions)),
                                   n25=length(intersect(D2_sig_positions, D5_sig_positions)),
                                   n34=length(intersect(D3_sig_positions, D4_sig_positions)),
                                   n35=length(intersect(D3_sig_positions, D5_sig_positions)),
                                   n45=length(intersect(D4_sig_positions, D5_sig_positions)),
                                   n123=length(intersect(intersect(D1_sig_positions, D2_sig_positions), D3_sig_positions)),
                                   n124=length(intersect(intersect(D1_sig_positions, D2_sig_positions), D4_sig_positions)),
                                   n125=length(intersect(intersect(D1_sig_positions, D2_sig_positions), D5_sig_positions)),
                                   n134=length(intersect(intersect(D1_sig_positions, D3_sig_positions), D4_sig_positions)),
                                   n135=length(intersect(intersect(D1_sig_positions, D3_sig_positions), D5_sig_positions)),
                                   n145=length(intersect(intersect(D1_sig_positions, D4_sig_positions), D5_sig_positions)),
                                   n234=length(intersect(intersect(D2_sig_positions, D3_sig_positions), D4_sig_positions)),
                                   n235=length(intersect(intersect(D2_sig_positions, D3_sig_positions), D5_sig_positions)),
                                   n245=length(intersect(intersect(D2_sig_positions, D4_sig_positions), D5_sig_positions)),
                                   n345=length(intersect(intersect(D3_sig_positions, D4_sig_positions), D5_sig_positions)),
                                   n1234=length(intersect(intersect(D1_sig_positions, D2_sig_positions), 
                                                          intersect(D3_sig_positions, D4_sig_positions))),
                                   n1235=length(intersect(intersect(D1_sig_positions, D2_sig_positions),
                                                          intersect(D3_sig_positions, D5_sig_positions))),
                                   n1245=length(intersect(intersect(D1_sig_positions, D2_sig_positions),
                                                          intersect(D4_sig_positions, D5_sig_positions))),
                                   n1345=length(intersect(intersect(D1_sig_positions, D3_sig_positions),
                                                          intersect(D4_sig_positions, D5_sig_positions))),
                                   n2345=length(intersect(intersect(D2_sig_positions, D3_sig_positions),
                                                          intersect(D4_sig_positions, D5_sig_positions))),
                                   n12345=length(intersect(intersect(D1_sig_positions, D2_sig_positions), 
                                                           intersect(intersect(D3_sig_positions, D4_sig_positions), D5_sig_positions))),
                                   category=c("D1", "D2", "D3", "D4", "D5"))

#############################
# Gene Venn Diagram objects #
#############################


C_venn_genes <- draw.quintuple.venn(area1=length(C1_sig_gene_names), area2=length(C2_sig_gene_names), area3=length(C3_sig_gene_names),
                    area4=length(C4_sig_gene_names), area5=length(C5_sig_gene_names),
                    n12=length(intersect(C1_sig_gene_names, C2_sig_gene_names)), 
                    n13=length(intersect(C1_sig_gene_names, C3_sig_gene_names)),
                    n14=length(intersect(C1_sig_gene_names, C4_sig_gene_names)),
                    n15=length(intersect(C1_sig_gene_names, C5_sig_gene_names)),
                    n23=length(intersect(C2_sig_gene_names, C3_sig_gene_names)),
                    n24=length(intersect(C2_sig_gene_names, C4_sig_gene_names)),
                    n25=length(intersect(C2_sig_gene_names, C5_sig_gene_names)),
                    n34=length(intersect(C3_sig_gene_names, C4_sig_gene_names)),
                    n35=length(intersect(C3_sig_gene_names, C5_sig_gene_names)),
                    n45=length(intersect(C4_sig_gene_names, C5_sig_gene_names)),
                    n123=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names), C3_sig_gene_names)),
                    n124=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names), C4_sig_gene_names)),
                    n125=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names), C5_sig_gene_names)),
                    n134=length(intersect(intersect(C1_sig_gene_names, C3_sig_gene_names), C4_sig_gene_names)),
                    n135=length(intersect(intersect(C1_sig_gene_names, C3_sig_gene_names), C5_sig_gene_names)),
                    n145=length(intersect(intersect(C1_sig_gene_names, C4_sig_gene_names), C5_sig_gene_names)),
                    n234=length(intersect(intersect(C2_sig_gene_names, C3_sig_gene_names), C4_sig_gene_names)),
                    n235=length(intersect(intersect(C2_sig_gene_names, C3_sig_gene_names), C5_sig_gene_names)),
                    n245=length(intersect(intersect(C2_sig_gene_names, C4_sig_gene_names), C5_sig_gene_names)),
                    n345=length(intersect(intersect(C3_sig_gene_names, C4_sig_gene_names), C5_sig_gene_names)),
                    n1234=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names), 
                                           intersect(C3_sig_gene_names, C4_sig_gene_names))),
                    n1235=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names),
                                           intersect(C3_sig_gene_names, C5_sig_gene_names))),
                    n1245=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names),
                                           intersect(C4_sig_gene_names, C5_sig_gene_names))),
                    n1345=length(intersect(intersect(C1_sig_gene_names, C3_sig_gene_names),
                                           intersect(C4_sig_gene_names, C5_sig_gene_names))),
                    n2345=length(intersect(intersect(C2_sig_gene_names, C3_sig_gene_names),
                                           intersect(C4_sig_gene_names, C5_sig_gene_names))),
                    n12345=length(intersect(intersect(C1_sig_gene_names, C2_sig_gene_names), 
                                            intersect(intersect(C3_sig_gene_names, C4_sig_gene_names), C5_sig_gene_names))),
                    category=c("C1", "C2", "C3", "C4", "C5"))





 
D_venn_genes <- draw.quintuple.venn(area1=length(D1_sig_gene_names), area2=length(D2_sig_gene_names), area3=length(D3_sig_gene_names),
                                    area4=length(D4_sig_gene_names), area5=length(D5_sig_gene_names),
                                    n12=length(intersect(D1_sig_gene_names, D2_sig_gene_names)), 
                                    n13=length(intersect(D1_sig_gene_names, D3_sig_gene_names)),
                                    n14=length(intersect(D1_sig_gene_names, D4_sig_gene_names)),
                                    n15=length(intersect(D1_sig_gene_names, D5_sig_gene_names)),
                                    n23=length(intersect(D2_sig_gene_names, D3_sig_gene_names)),
                                    n24=length(intersect(D2_sig_gene_names, D4_sig_gene_names)),
                                    n25=length(intersect(D2_sig_gene_names, D5_sig_gene_names)),
                                    n34=length(intersect(D3_sig_gene_names, D4_sig_gene_names)),
                                    n35=length(intersect(D3_sig_gene_names, D5_sig_gene_names)),
                                    n45=length(intersect(D4_sig_gene_names, D5_sig_gene_names)),
                                    n123=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names), D3_sig_gene_names)),
                                    n124=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names), D4_sig_gene_names)),
                                    n125=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names), D5_sig_gene_names)),
                                    n134=length(intersect(intersect(D1_sig_gene_names, D3_sig_gene_names), D4_sig_gene_names)),
                                    n135=length(intersect(intersect(D1_sig_gene_names, D3_sig_gene_names), D5_sig_gene_names)),
                                    n145=length(intersect(intersect(D1_sig_gene_names, D4_sig_gene_names), D5_sig_gene_names)),
                                    n234=length(intersect(intersect(D2_sig_gene_names, D3_sig_gene_names), D4_sig_gene_names)),
                                    n235=length(intersect(intersect(D2_sig_gene_names, D3_sig_gene_names), D5_sig_gene_names)),
                                    n245=length(intersect(intersect(D2_sig_gene_names, D4_sig_gene_names), D5_sig_gene_names)),
                                    n345=length(intersect(intersect(D3_sig_gene_names, D4_sig_gene_names), D5_sig_gene_names)),
                                    n1234=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names), 
                                                           intersect(D3_sig_gene_names, D4_sig_gene_names))),
                                    n1235=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names),
                                                           intersect(D3_sig_gene_names, D5_sig_gene_names))),
                                    n1245=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names),
                                                           intersect(D4_sig_gene_names, D5_sig_gene_names))),
                                    n1345=length(intersect(intersect(D1_sig_gene_names, D3_sig_gene_names),
                                                           intersect(D4_sig_gene_names, D5_sig_gene_names))),
                                    n2345=length(intersect(intersect(D2_sig_gene_names, D3_sig_gene_names),
                                                           intersect(D4_sig_gene_names, D5_sig_gene_names))),
                                    n12345=length(intersect(intersect(D1_sig_gene_names, D2_sig_gene_names), 
                                                            intersect(intersect(D3_sig_gene_names, D4_sig_gene_names), D5_sig_gene_names))),
                                    category=c("D1", "D2", "D3", "D4", "D5"))



#####################
# Save images on a  #
# single page for   #
# Figure 1          #
#####################

setwd('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/gene_list_overlap_testing')

pdf(file="Venn_SNP_and_gene_overlap.pdf", width=12, height=12)
grid.layout(nrow=2, ncol=2)
grid.draw(C_venn_SNPs)
grid.text(label="SNPs", x=0.9, y=0.9)
grid.newpage()
grid.draw(C_venn_genes)
grid.text(label="Genes", x=0.9, y=0.9)
grid.newpage()
grid.draw(D_venn_SNPs)
grid.text(label="SNPs", x=0.9, y=0.9)
grid.newpage()
grid.draw(D_venn_genes)
grid.text(label="Genes", x=0.9, y=0.9)
dev.off()

##########################
# Heatmap visualisations #
##########################

C_gene_names <- union(union(union(C1_sig_gene_names, 
                                  C2_sig_gene_names), 
                            C3_sig_gene_names), 
                      union(C4_sig_gene_names, 
                            C5_sig_gene_names))


C_genes_in_D1 <- intersect(D1_withC_gene_names, C_gene_names)
C_genes_in_D2 <- intersect(D2_withC_gene_names, C_gene_names)
C_genes_in_D3 <- intersect(D3_withC_gene_names, C_gene_names)
C_genes_in_D4 <- intersect(D4_withC_gene_names, C_gene_names)
C_genes_in_D5 <- intersect(D5_withC_gene_names, C_gene_names)

C_gene_pres_abs <- data.frame(row.names=C_gene_names,
                              C1=C_gene_names%in%C1_sig_gene_names,
                              C2=C_gene_names%in%C2_sig_gene_names,
                              C3=C_gene_names%in%C3_sig_gene_names,
                              C4=C_gene_names%in%C4_sig_gene_names,
                              C5=C_gene_names%in%C5_sig_gene_names,
                              D1=C_gene_names%in%D1_withC_gene_names,
                              D2=C_gene_names%in%D2_withC_gene_names,
                              D3=C_gene_names%in%D3_withC_gene_names,
                              D4=C_gene_names%in%D4_withC_gene_names,
                              D5=C_gene_names%in%D5_withC_gene_names)
C_gene_pres_abs_sorted <- C_gene_pres_abs[order(ordered(row.names(C_gene_pres_abs))),]

# Obtain the genomic location of the genes in the list for sorting
# 
# gene_list_for_conversion <- read.table("/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/All_sig_gene_list_for_name_conversion.txt",
#                                        header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="\"")
# 


C_gene_pres_abs_sorted <- t(as.matrix(C_gene_pres_abs_sorted, row.names=row.names(C_gene_pres_abs_sorted), col.names=c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5")))

pdf("C gene overlap heatmap visualisation 150728.pdf", height=50, width=5)
par(mar = c(5,5,4,2) + 0.1, mgp=c(3,0.2,0))
image(z = C_gene_pres_abs_sorted, col = c("grey", "red"), axes = FALSE)
axis(side = 2, labels = colnames(C_gene_pres_abs_sorted), 
     las=2, at=seq(0, 1, by=1/(ncol(C_gene_pres_abs_sorted)-1)), 
     cex.axis=0.15, tick=FALSE)
axis(side = 1, labels = rownames(C_gene_pres_abs_sorted), 
     at = seq(0, 1, by=1/(nrow(C_gene_pres_abs_sorted)-1)), 
     las = 1, tick=FALSE)
dev.off()


D_gene_names <- union(union(union(D1_sig_gene_names, 
                                  D2_sig_gene_names), 
                            D3_sig_gene_names), 
                      union(D4_sig_gene_names, 
                            D5_sig_gene_names))


D_gene_pres_abs <- data.frame(row.names=D_gene_names,
                              D1=D_gene_names%in%D1_sig_gene_names,
                              D2=D_gene_names%in%D2_sig_gene_names,
                              D3=D_gene_names%in%D3_sig_gene_names,
                              D4=D_gene_names%in%D4_sig_gene_names,
                              D5=D_gene_names%in%D5_sig_gene_names)
D_gene_pres_abs_sorted <- D_gene_pres_abs[order(ordered(row.names(D_gene_pres_abs))),]


D_gene_pres_abs_sorted <- t(as.matrix(D_gene_pres_abs_sorted, 
                                      row.names=row.names(D_gene_pres_abs_sorted), 
                                      col.names=c("D1", "D2", "D3", "D4", "D5")))



pdf("D gene overlap without C genes heatmap visualisation 150728.pdf", height=120, width=2)
par(mar = c(5,5,4,2) + 0.1, mgp=c(3,0.2,0))
image(z = D_gene_pres_abs_sorted, col = c("grey", "green"), axes = FALSE)
axis(side = 2, labels = colnames(D_gene_pres_abs_sorted), 
     las=2, at=seq(0, 1, by=1/(ncol(D_gene_pres_abs_sorted)-1)), 
     cex.axis=0.05, tick=FALSE)
axis(side = 1, labels = rownames(D_gene_pres_abs_sorted), 
     at = seq(0, 1, by=1/(nrow(D_gene_pres_abs_sorted)-1)), 
     las = 1, tick=FALSE, cex.axis=0.25)
dev.off()


########
########
########



# D_venn_GO <- draw.quintuple.venn(area1=length(D1_sig_GO), area2=length(D2_sig_GO), area3=length(D3_sig_GO),
#                                  area4=length(D4_sig_GO), area5=length(D5_sig_GO),
#                                  n12=length(intersect(D1_sig_GO, D2_sig_GO)), 
#                                  n13=length(intersect(D1_sig_GO, D3_sig_GO)),
#                                  n14=length(intersect(D1_sig_GO, D4_sig_GO)),
#                                  n15=length(intersect(D1_sig_GO, D5_sig_GO)),
#                                  n23=length(intersect(D2_sig_GO, D3_sig_GO)),
#                                  n24=length(intersect(D2_sig_GO, D4_sig_GO)),
#                                  n25=length(intersect(D2_sig_GO, D5_sig_GO)),
#                                  n34=length(intersect(D3_sig_GO, D4_sig_GO)),
#                                  n35=length(intersect(D3_sig_GO, D5_sig_GO)),
#                                  n45=length(intersect(D4_sig_GO, D5_sig_GO)),
#                                  n123=length(intersect(intersect(D1_sig_GO, D2_sig_GO), D3_sig_GO)),
#                                  n124=length(intersect(intersect(D1_sig_GO, D2_sig_GO), D4_sig_GO)),
#                                  n125=length(intersect(intersect(D1_sig_GO, D2_sig_GO), D5_sig_GO)),
#                                  n134=length(intersect(intersect(D1_sig_GO, D3_sig_GO), D4_sig_GO)),
#                                  n135=length(intersect(intersect(D1_sig_GO, D3_sig_GO), D5_sig_GO)),
#                                  n145=length(intersect(intersect(D1_sig_GO, D4_sig_GO), D5_sig_GO)),
#                                  n234=length(intersect(intersect(D2_sig_GO, D3_sig_GO), D4_sig_GO)),
#                                  n235=length(intersect(intersect(D2_sig_GO, D3_sig_GO), D5_sig_GO)),
#                                  n245=length(intersect(intersect(D2_sig_GO, D4_sig_GO), D5_sig_GO)),
#                                  n345=length(intersect(intersect(D3_sig_GO, D4_sig_GO), D5_sig_GO)),
#                                  n1234=length(intersect(intersect(D1_sig_GO, D2_sig_GO), 
#                                                         intersect(D3_sig_GO, D4_sig_GO))),
#                                  n1235=length(intersect(intersect(D1_sig_GO, D2_sig_GO),
#                                                         intersect(D3_sig_GO, D5_sig_GO))),
#                                  n1245=length(intersect(intersect(D1_sig_GO, D2_sig_GO),
#                                                         intersect(D4_sig_GO, D5_sig_GO))),
#                                  n1345=length(intersect(intersect(D1_sig_GO, D3_sig_GO),
#                                                         intersect(D4_sig_GO, D5_sig_GO))),
#                                  n2345=length(intersect(intersect(D2_sig_GO, D3_sig_GO),
#                                                         intersect(D4_sig_GO, D5_sig_GO))),
#                                  n12345=length(intersect(intersect(D1_sig_GO, D2_sig_GO), 
#                                                          intersect(intersect(D3_sig_GO, D4_sig_GO), D5_sig_GO))),
#                                  category=c("D1", "D2", "D3", "D4", "D5"))
# 
# pdf(file="Overlap_for_D_replicates_GO_categories.pdf", width=6, height=6)
# par(mfrow=c(1,2))
# grid.draw(D_venn_GO)
# grid.text(label="GO Categories", x=0.9, y=0.9)
# dev.off()
# 
# 
# D_GO_names <- union(union(union(D1_sig_GO, D2_sig_GO), D3_sig_GO), union(D4_sig_GO, D5_sig_GO))
# 
# D_GO_pres_abs <- data.frame(row.names=D_GO_names,
#                             D1=D_GO_names%in%D1_sig_GO,
#                             D2=D_GO_names%in%D2_sig_GO,
#                             D3=D_GO_names%in%D3_sig_GO,
#                             D4=D_GO_names%in%D4_sig_GO,
#                             D5=D_GO_names%in%D5_sig_GO)
# D_GO_pres_abs <- t(as.matrix(D_GO_pres_abs, row.names=D_GO_names, col.names=c("D1", "D2", "D3", "D4", "D5")))
# 
# 
# pdf("D GO overlap without C genes heatmap visualisation.pdf", height=50, width=6)
# par(mar = c(5,9,4,2) + 0.1, mgp=c(3,0.2,0))
# image(z = D_GO_pres_abs, col = c("grey", "blue"), axes = FALSE)
# axis(side = 2, labels = colnames(D_GO_pres_abs), las=2, at=seq(0, 1, by=1/(ncol(D_GO_pres_abs)-1)), cex.axis=0.3, tick=FALSE)
# axis(side = 1, labels = rownames(D_GO_pres_abs), at = seq(0, 1, by=1/(nrow(D_GO_pres_abs)-1)), las = 1, tick=FALSE)
# dev.off()

