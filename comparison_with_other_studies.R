obs <- readRDS(file="/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/gene_list_overlap_testing/Gene_presence_count_in_simulated_gene_lists_150826.rds")

get_D_sim_counts <- function(gene_name, count_table){
  temp_values <- count_table[rownames(count_table)==gene_name, 6:10]
  return(temp_values/1000)
}

get_D_presence_vector <- function(gene_name, list_of_gene_lists){
  presence_vector <- c()
  for(i in list_of_gene_lists){
    temp_presence <- gene_name%in%i
    presence_vector <- c(presence_vector, temp_presence)
  }
  names(presence_vector) <- names(list_of_gene_lists)
  return(presence_vector)
}

gwas_and_mbe_genes <- c('dnc', 'ct', 'Rbf', 'sdk', 'Pde9', 'CG43658', 
                        'CrebB', 'santa-maria', 'Nckx30C', 'Trim9', 'nub', 
                        'sick', 'mam', 'CG8910', 'Amph', 'CG1688', 'jbug', 
                        'fab1', 'CG9864', 'luna', 'Fili', 'Strn-Mlck', 'fz', 
                        'nkd', 'klu', 'CG11486', 'CG33275', 'sowah', 'CG18769', 
                        'CG32264', 'Rbp6', 'Ect4', 'Abd-B', 'slo', 'fru', 'klg', 
                        'GluClalpha', 'CG4390', 'CG7084', 'Mtl', 'CG13643', 
                        'CG31176', 'cv-c', 'Gprk2', 'SNF4Agamma')

output_matrix <- matrix(NA, nrow=length(gwas_and_mbe_genes), ncol=5, byrow=TRUE)
rownames(output_matrix) <- gwas_and_mbe_genes

for(i in 1:length(gwas_and_mbe_genes)){
  temp_gene <- gwas_and_mbe_genes[i]
  temp_counts <- get_D_sim_counts(gene_name=temp_gene, count_table=obs)
  if(length(temp_counts)>0){
    output_matrix[i,] <- temp_counts
  }
}

lit_genes <- c('Trpm', 'dnc', 'Rel', 'msn', 'aop', 'nompC', 'CG33958', 
               'shark', 'Pkd2', 'Pde1c', 'Src64B', 'Rho1', 'trp', 'Nplp4', 
               'Nos', 'Dok', 'Ilp1', 'Aplip1', 'Tak1', 'Jra', 'Desat1', 
               'Trpml', 'dock', 'nan', 'wtrw', 'Btk29A', 'CapaR', 'Pak', 
               'iav', 'CG34357', 'hep', 'Atf-2', 'CG5171', 'CG3523', 'bsk', 
               'Trpgamma', 'Rac2', 'Nplp3', 'InR', 'pyd', 'egr', 'Ilp4', 
               'pain', 'Fatp', 'Gyc76C', 'Capa', 'Ilp3', 'TrpA1', 'Pgm', 
               'trpl', 'ITP', 'Src42A', 'kay', 'slpr', 'medo', 'Tk', 'Cdc42', 
               'for', 'puc', 'Desat2', 'Nplp2', 'CG6262', 'CYLD', 'cno', 
               'Ilp2', 'Gadd45', 'Mtl', 'sNPF', 'Pde9', 'chico', 'Cka', 
               'CG5177', 'Tps1', 'Nplp1', 'Ilp6', 'Pde6', 'hppy', 'poly', 
               'klu', 'v_2_k05816', 'Ilp8', 'pyx', 'crc', 'cbt', 'Treh', 'Cyp4g1', 
               'Dh44-R2', 'Ilp7', 'Mkk4', 'Ilp5', 'Pde11', 'Rac1', 'Pkg21D')



setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

# Obtain full list of sig C genes

C_sig_genes <- read.csv('Putative_lab_adaptation_gene_list.txt', header=FALSE,
                        stringsAsFactors=FALSE)[,1]

# Now get D1, D4, D5 gene lists

D1 <- read.table(file="D1_noC_sig_gene_list.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1]
D2 <- read.table(file="D2_noC_sig_gene_list.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1]
D3 <- read.table(file="D3_noC_sig_gene_list.txt", sep="\t", 
           header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1]
D4 <- read.table(file="D4_noC_sig_gene_list.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1]
D5 <- read.table(file="D5_noC_sig_gene_list.txt", sep="\t", 
                 header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1]

lists_to_check <- list(D1=D1, D2=D2, D3=D3, D4=D4, D5=D5, allC=C_sig_genes)

get_D_presence_vector('for', lists_to_check)

output_matrix_2 <- matrix(NA, nrow=length(lit_genes), ncol=11, byrow=TRUE)
rownames(output_matrix_2) <- lit_genes

for(i in 1:length(lit_genes)){
  temp_gene <- lit_genes[i]
  temp_presence <- get_D_presence_vector(temp_gene, lists_to_check)
  temp_counts <- get_D_sim_counts(gene_name=temp_gene, count_table=obs)
  output_matrix_2[i,1:6] <- temp_presence
  if(length(temp_counts)>0){
    output_matrix_2[i,7:11] <- temp_counts
  }
}

#import background network...

library(igraph)

#############
# FUNCTIONS #
#############

make_zero_degree_subgraph <- function(background_network, gene_list){
  gene_present <- V(background_network)$name%in%gene_list
  numbering <- 1:length(V(background_network)$name)
  numbers_to_keep <- numbering[gene_present]
  pre_zero_subgraph <- induced.subgraph(background_network, numbers_to_keep)
  nonzero_degree_vertices <- degree(pre_zero_subgraph)>0
  numbering_2 <- 1:length(V(pre_zero_subgraph)$name)
  numbers_to_keep_2 <- numbering_2[nonzero_degree_vertices]
  zero_subgraph <- induced.subgraph(pre_zero_subgraph, numbers_to_keep_2)
  zero_subgraph <- simplify(zero_subgraph, remove.multiple=TRUE, remove.loops=TRUE)
  edgelist <- get.edgelist(zero_subgraph)
  edgelist2 <- data.frame(paste(edgelist[,1], edgelist[,2], sep="_"), paste(edgelist[,2], edgelist[,1], sep="_"))
  zero_subgraph <- set.edge.attribute(zero_subgraph, "name1", value=as.character(edgelist2[,1]))
  zero_subgraph <- set.edge.attribute(zero_subgraph, "name2", value=as.character(edgelist2[,2]))
  
  return(zero_subgraph)
}

make_first_degree_subgraph <- function(background_network, gene_list){
  background_edge_list <- get.edgelist(background_network)
  #gene_present <- which(background_edge_list[,1]%in%gene_list|background_edge_list[,2]%in%gene_list)
  #pre_first_subgraph <- subgraph.edges(background_network, eids=gene_present)
  pre_first_subgraph <- induced.subgraph(graph=background_network, 
                                         vids=unlist(neighborhood(graph=background_network,order=1, 
                                                                  nodes=which(V(background_network)$name%in%gene_list))))
#   nonzero_degree_vertices <- degree(pre_first_subgraph)>0
#   numbering_2 <- 1:length(V(pre_first_subgraph)$name)
#   numbers_to_keep_2 <- numbering_2[nonzero_degree_vertices]
#   first_subgraph <- induced.subgraph(pre_first_subgraph, numbers_to_keep_2)
  first_subgraph <- simplify(pre_first_subgraph, remove.multiple=TRUE, remove.loops=TRUE)
  edgelist <- get.edgelist(first_subgraph)
  edgelist2 <- data.frame(paste(edgelist[,1], edgelist[,2], sep="_"), paste(edgelist[,2], edgelist[,1], sep="_"))
  first_subgraph <- set.edge.attribute(first_subgraph, "name1", value=as.character(edgelist2[,1]))
  first_subgraph <- set.edge.attribute(first_subgraph, "name2", value=as.character(edgelist2[,2]))
  
  return(first_subgraph)
}

make_intersection <- function(graph1, graph2){
  graph1_graph2_intersection <- graph.intersection(graph1, graph2)
  graph1_graph2 <- delete.vertices(graph1_graph2_intersection, 
                                   which(degree(graph1_graph2_intersection)<1))
  return(graph1_graph2)
}





#First step: to get a good background Drosophila protein-protein interaction network.

#Avoiding interologs (interaction info from other species) on Melissa Davis' 
#recommendation.

#using all info from DroID as of 29th May 2015 (v 2014_10)

setwd("~/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment/DroID_v2014_10")

P1_list <- c()
P2_list <- c()
for(i in list.files()){
  if (i!="not_using"){
    temp_list <- read.csv(i, sep="\t", stringsAsFactors=FALSE)
    P1_list <- c(P1_list, temp_list[,1])
    P2_list <- c(P2_list, temp_list[,2])
  }
}
int_list <- data.frame(P1_list, P2_list)

int_list <- int_list[duplicated(int_list)==FALSE,]


test_network <- graph.data.frame(int_list, directed=FALSE, vertices=NULL)
test_network <- simplify(test_network, remove.multiple=TRUE, remove.loops=TRUE)
test_network <- delete.vertices(test_network, which(degree(test_network)==0))

#can't access edge names!!! gah
#have to interpret somehow from E(test_network)

edgelist <- get.edgelist(test_network)
edgelist2 <- data.frame(paste(edgelist[,1], edgelist[,2], sep="_"), paste(edgelist[,2], edgelist[,1], sep="_"))
colnames(edgelist2) <- c("order1", "order2")

test_network <- set.edge.attribute(test_network, "name1", value=as.character(edgelist2[,1]))
test_network <- set.edge.attribute(test_network, "name2", value=as.character(edgelist2[,2]))


vcount(test_network)
ecount(test_network)

# get table of presence and absence of physiology candidates in gene lists, along with FBgn names

phys <- read.table(file="/Users/pgriffin/Documents/Drosophila Selection Experiment/Investigating_physiology_candidates.txt",
                   stringsAsFactors=FALSE, header=TRUE, sep="\t")


levels(as.factor(phys$Category_1))

phys_FBgn <- phys$FBgn

phys_vertices <- which(V(test_network)$name%in%phys_FBgn)

phys_network <- induced.subgraph(test_network, v=phys_vertices)
phys_1deg_network <- make_first_degree_subgraph(test_network, phys_FBgn)

plot(phys_network, vertex.label.cex=0.2)



#colour vertex lines by category
cols <- c("blue", "forestgreen", "black", "purple", "orange", "yellow", "grey", "skyblue")
line_col <- cols[as.factor(phys_sorted$Category_1)]

#make a good basic layout
#layout1 <- layout.fruchterman.reingold(phys_network)
layout1 <- readRDS(file="/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/physiology_candidate_network_layout.rds")

unconnected_vertex_indices <- which(degree(phys_network)==0)
cat1_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="Calcium signalling")
cat2_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="cGMP or cAMP signalling")
cat3_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="Cuticular hydrocarbons")
cat4_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="Desiccation specific")
cat5_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="Insulin receptor binding")
cat6_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="MAPK or JNK")
cat7_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="Trehalose metabolism")
cat8_to_change <- which(degree(phys_network)==0&phys_sorted[,2]=="TRP ion channels")


layout1[cat1_to_change,2] <- -60
layout1[cat1_to_change,1] <- 18

layout1[cat2_to_change,2] <- c(rep(seq(-36, 0, by=12), times=3), -36)
layout1[cat2_to_change,1] <- c(rep(seq(30, 54, by=12), times=4), 66)
layout1[cat3_to_change,2] <- c(-72, -60, -72, -60, -48)
layout1[cat3_to_change,1] <- c(-6, 6, 6, -6, -6)
layout1[cat4_to_change,2] <- c(12, 24, 24)
layout1[cat4_to_change,1] <- c(30, 30, 18)
layout1[cat5_to_change,2] <- c(36, 48, 36, 48)
layout1[cat5_to_change,1] <- c(6, 6, 18, 18)
layout1[cat6_to_change,2] <- c(-36, -36, -36, -48, -48)
layout1[cat6_to_change,1] <- c(-6, 6, 18, 6, 18)
layout1[cat7_to_change,2] <- c(-48, -48, -48, -60, -60, -60)
layout1[cat7_to_change,1] <- c(30, 42, 54, 30, 42, 54)
layout1[cat8_to_change,2] <- c(36, 48, 48, 48)
layout1[cat8_to_change,1] <- c(30, 30, 42, 54)

saveRDS(layout1, file="/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/physiology_candidate_network_layout.rds")

#create a new vertex attribute to reflect the presence of each gene in each
#dataset. 
# This 'presence' value is 1 - the proportion of resampling iterations where the gene was detected
# (multiplied by 0 if the gene was not a hit in the observed data)
# so basically genes that are present and 'less likely to be a chance result' are red,
# genes that are present, but likely to be a chance result, are a lighter shade
# genes that are absent are white

phys_values <- phys[,c("FBgn", "Category_1", "d1_scale", "d2_scaled", "d3_scaled", "d4_scaled", "d5_scaled", "C_presence")]

phys_sorted <- phys_values[match(get.vertex.attribute(phys_network, "name"), phys_values[,1]),]


D1_colour <- c()
for(i in phys_sorted[,'d1_scale']){
  temp_colour <- rgb(1, 0, 0, i, maxColorValue=1)
  D1_colour <- c(D1_colour, temp_colour)
}

D2_colour <- c()
for(i in phys_sorted[,'d2_scaled']){
  temp_colour <- rgb(1, 0, 0, i, maxColorValue=1)
  D2_colour <- c(D2_colour, temp_colour)
}

D3_colour <- c()
for(i in phys_sorted[,'d3_scaled']){
  temp_colour <- rgb(1, 0, 0, i, maxColorValue=1)
  D3_colour <- c(D3_colour, temp_colour)
}

D4_colour <- c()
for(i in phys_sorted[,'d4_scaled']){
  temp_colour <- rgb(1, 0, 0, i, maxColorValue=1)
  D4_colour <- c(D4_colour, temp_colour)
}

D5_colour <- c()
for(i in phys_sorted[,'d5_scaled']){
  temp_colour <- rgb(1, 0, 0, i, maxColorValue=1)
  D5_colour <- c(D5_colour, temp_colour)
}

D1_colour[which(phys_sorted$C_presence==1)] <- rgb(0, 0, 0, 0.5, maxColorValue=1)
D2_colour[which(phys_sorted$C_presence==1)] <- rgb(0, 0, 0, 0.5, maxColorValue=1)
D3_colour[which(phys_sorted$C_presence==1)] <- rgb(0, 0, 0, 0.5, maxColorValue=1)
D4_colour[which(phys_sorted$C_presence==1)] <- rgb(0, 0, 0, 0.5, maxColorValue=1)
D5_colour[which(phys_sorted$C_presence==1)] <- rgb(0, 0, 0, 0.5, maxColorValue=1)


pdf(file="Physiological Candidate Representation.pdf", width=28, height=28)
par(mfrow=c(3, 2))
for(i in c('D1', 'D2', 'D3', 'D4', 'D5')){
  temp_colour <- get(paste(i, "_colour", sep=""))
  V(phys_network)$color <- temp_colour
  plot(phys_network, layout=layout1, vertex.label.cex=0.2, 
       vertex.size=8, vertex.frame.color=line_col)
}
dev.off()


gwas_phys <- c('dnc', 'Pde9', 'klu', 'aop', 'Dok', 
               'hep', 'kay', 'msn', 'Mtl', 'puc', 'Rac1', 
               'Src64B', 'InR', 'CG3523', 'Treh', 'Pgm')
gwas_FBgn <- phys[match(gwas_phys, phys$name),"FBgn"]
mbe_phys <- c('trpl', 'CG34357', 'dnc', 'Pde1c', 'Pde9', 
              'Pde11', 'Nos', 'klu', 'ITP', 'sNPF', 'Mtl', 
              'slpr', 'poly')
mbe_FBgn <- phys[match(mbe_phys, phys$name),"FBgn"]

gwas_vertices <- which(V(phys_network)$name%in%gwas_FBgn)

mbe_vertices <- which(V(phys_network)$name%in%mbe_FBgn)

pdf(file="Physiological Candidate Representation for GWAS and MBE.pdf", width=28, height=28)
par(mfrow=c(1, 2))
for(i in list(gwas_vertices, mbe_vertices)){
  #temp_colour <- get(paste(i, "_colour", sep=""))
  V(phys_network)$color <- "white"
  V(phys_network)$color[i] <- "red"
  plot(phys_network, layout=layout1, vertex.label.cex=0.2, 
       vertex.size=8, vertex.frame.color=line_col)
}
dev.off()

#leave first level where it is

plot(phys_network, layout=layout1, vertex.label.cex=0.2, 
     vertex.size=8, vertex.frame.color=line_col)







cuticular <- subset(phys, phys$Category_1=="Cuticular hydrocarbons" | phys$Category_2=="Cuticular hydrocarbons")

cuticular_FBgn <- cuticular$FBgn

cuticular_vertices <- which(V(test_network)$name%in%cuticular_FBgn)

cuticular_network <- induced.subgraph(test_network, v=cuticular_vertices)


mapkjnk <- subset(phys, phys$Category_1=="MAPK or JNK" | phys$Category_2=="MAPK or JNK")

mapkjnk_FBgn <- mapkjnk$FBgn

mapkjnk_vertices <- which(V(test_network)$name%in%mapkjnk_FBgn)

mapkjnk_network <- induced.subgraph(test_network, v=mapkjnk_vertices)
