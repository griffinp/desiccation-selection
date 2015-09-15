
library(VennDiagram)
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

make_intersection <- function(graph1, graph2){
  graph1_graph2_intersection <- graph.intersection(graph1, graph2)
  graph1_graph2 <- delete.vertices(graph1_graph2_intersection, 
                                   which(degree(graph1_graph2_intersection)<1))
  return(graph1_graph2)
}

########
# MAIN #
########



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

#Okay, this is a good background network.

# Alternative background network: the one 

# setwd("~/Documents/Drosophila Selection Experiment/network_building_and_overlap_testing")
# rahul_network <- read.graph(file = "Flynet_summary_digraph.graphml", format = "graphml")
# rahul_network <- simplify(rahul_network, remove.multiple=TRUE, remove.loops=TRUE)
# rahul_network <- delete.vertices(rahul_network, which(degree(rahul_network)==0))
# #r_edgelist <- get.edgelist(rahul_network)
# r_names <- get.vertex.attribute(rahul_network, "id")
# V(rahul_network)$name <- r_names
# #set.vertex.attribute(rahul_network, "name", r_names)
# 
# test_network <- rahul_network
# 

#Things to do next:
#1) import a background gene list (i.e. all genes in Drosophila genome)
#2) convert my gene lists into FBgn format to match properly
#with the network data (doing this through Flybase with manual checking)
#3) clearly set up hypotheses for testing and figure out what background
#is appropriate for resampling from for each
#4) to make subgraphs of the background network, use the V(g)$name access
#approach (vertices are otherwise identified by numbers, which can change)
setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")

D1_list <- read.csv("D1_noC_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
D1_list <- unique(D1_list[,1])
D2_list <- read.csv("D2_noC_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
D2_list <- unique(D2_list[,1])
D3_list <- read.csv("D3_noC_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
D3_list <- unique(D3_list[,1])
D4_list <- read.csv("D4_noC_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
D4_list <- unique(D4_list[,1])
D5_list <- read.csv("D5_noC_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
D5_list <- unique(D5_list[,1])

C1_list <- read.csv("C1_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
C1_list <- unique(C1_list[,1])
C2_list <- read.csv("C2_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
C2_list <- unique(C2_list[,1])
C3_list <- read.csv("C3_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
C3_list <- unique(C3_list[,1])
C4_list <- read.csv("C4_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
C4_list <- unique(C4_list[,1])
C5_list <- read.csv("C5_sig_gene_names_FBgn_format.txt", header=FALSE, stringsAsFactors=FALSE)
C5_list <- unique(C5_list[,1])


D1_zero <- make_zero_degree_subgraph(test_network, D1_list)
D2_zero <- make_zero_degree_subgraph(test_network, D2_list)
D3_zero <- make_zero_degree_subgraph(test_network, D3_list)
D4_zero <- make_zero_degree_subgraph(test_network, D4_list)
D5_zero <- make_zero_degree_subgraph(test_network, D5_list)
C1_zero <- make_zero_degree_subgraph(test_network, C1_list)
C2_zero <- make_zero_degree_subgraph(test_network, C2_list)
C3_zero <- make_zero_degree_subgraph(test_network, C3_list)
C4_zero <- make_zero_degree_subgraph(test_network, C4_list)
C5_zero <- make_zero_degree_subgraph(test_network, C5_list)

allD_list <- unique(c(D1_list, D2_list, D3_list, D4_list, D5_list))
D_union <- make_zero_degree_subgraph(test_network, allD_list)



D1_D2 <- make_intersection(D1_zero, D2_zero)
D1_D3 <- make_intersection(D1_zero, D3_zero)
D1_D4 <- make_intersection(D1_zero, D4_zero)
D1_D5 <- make_intersection(D1_zero, D5_zero)
D2_D3 <- make_intersection(D2_zero, D3_zero)
D2_D4 <- make_intersection(D2_zero, D4_zero)
D2_D5 <- make_intersection(D2_zero, D5_zero)
D3_D4 <- make_intersection(D3_zero, D4_zero)
D3_D5 <- make_intersection(D3_zero, D5_zero)
D4_D5 <- make_intersection(D4_zero, D5_zero)
D1_D2_D3 <- make_intersection(D1_D2, D3_zero)
D1_D2_D4 <- make_intersection(D1_D2, D4_zero)
D1_D2_D5 <- make_intersection(D1_D2, D5_zero)
D1_D3_D4 <- make_intersection(D1_D3, D4_zero)
D1_D3_D5 <- make_intersection(D1_D3, D5_zero)
D1_D4_D5 <- make_intersection(D1_D4, D5_zero)
D2_D3_D4 <- make_intersection(D2_D3, D4_zero)
D2_D3_D5 <- make_intersection(D2_D3, D5_zero)
D2_D4_D5 <- make_intersection(D2_D4, D5_zero)
D3_D4_D5 <- make_intersection(D3_D4, D5_zero)
D1_D2_D3_D4 <- make_intersection(D1_D2, D3_D4)
D1_D2_D3_D5 <- make_intersection(D1_D2, D3_D5)
D1_D2_D4_D5 <- make_intersection(D1_D2, D4_D5)
D1_D3_D4_D5 <- make_intersection(D1_D3, D4_D5)
D2_D3_D4_D5 <- make_intersection(D2_D3, D4_D5)
D1_D2_D3_D4_D5 <- make_intersection(D1_D2_D3, D4_D5)

# edges_in_2 <- which(E(test_network)$name1 %in% E(D1_D2)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D3)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D5)$name1_1 | 
#                       E(test_network)$name1 %in% E(D2_D3)$name1_1 |
#                       E(test_network)$name1 %in% E(D2_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D2_D5)$name1_1 | 
#                       E(test_network)$name1 %in% E(D3_D4)$name1_1 |
#                       E(test_network)$name1 %in% E(D3_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D4_D5)$name1_1)                      )
# edges_in_3 <- which(E(test_network)$name1 %in% E(D1_D2_D3)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D2_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D2_D5)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D3_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D1_D3_D5)$name1_1 |
#                       E(test_network)$name1 %in% E(D1_D4_D5)$name1_1 | 
#                       E(test_network)$name1 %in% E(D2_D3_D4)$name1_1 | 
#                       E(test_network)$name1 %in% E(D2_D3_D5)$name1_1 |
#                       E(test_network)$name1 %in% E(D2_D4_D5)$name1_1 | 
#                       E(test_network)$name1 %in% E(D3_D4_D5)$name1_1)
# edges_in_4 <- which(E(test_network)$name1 %in% E(D1_D2_D3_D4)$name1_1_1 | 
#                       E(test_network)$name1 %in% E(D1_D2_D3_D5)$name1_1_1 | 
#                       E(test_network)$name1 %in% E(D1_D2_D4_D5)$name1_1_1 | 
#                       E(test_network)$name1 %in% E(D1_D3_D4_D5)$name1_1_1 | 
#                       E(test_network)$name1 %in% E(D2_D3_D4_D5)$name1_1_1)
# edges_in_5 <- which(E(test_network)$name1 %in% E(D1_D2_D3_D4_D5)$name1_1_1)
# 
# edge_overlap_colour <- rep(rgb(100, 100, 100, 100, maxColorValue=255), length=ecount(test_network))
# edge_overlap_colour[edges_in_2] <- rgb(200,100,0,100, maxColorValue=255)
# edge_overlap_colour[edges_in_3] <- rgb(255,140,0,100, maxColorValue=255)
# edge_overlap_colour[edges_in_4] <- rgb(0,255,0,100, maxColorValue=255)
# edge_overlap_colour[edges_in_5] <- rgb(0,0,255,100, maxColorValue=255)
# test_network <- set.edge.attribute(test_network, name="colour", value=edge_overlap_colour)


library(VennDiagram)

D_venn_edges <- draw.quintuple.venn(area1=ecount(D1_zero), area2=ecount(D2_zero), area3=ecount(D3_zero),
                                   area4=ecount(D4_zero), area5=ecount(D5_zero),
                                   n12=ecount(D1_D2), n13=ecount(D1_D3),
                                   n14=ecount(D1_D4), n15=ecount(D1_D5),
                                   n23=ecount(D2_D3), n24=ecount(D2_D4),
                                   n25=ecount(D2_D5), n34=ecount(D3_D4),
                                   n35=ecount(D3_D5), n45=ecount(D4_D5),
                                   n123=ecount(D1_D2_D3), n124=ecount(D1_D2_D4),
                                   n125=ecount(D1_D2_D5), n134=ecount(D1_D3_D4),
                                   n135=ecount(D1_D3_D5), n145=ecount(D1_D4_D5),
                                   n234=ecount(D2_D3_D4), n235=ecount(D2_D3_D5),
                                   n245=ecount(D2_D4_D5), n345=ecount(D3_D4_D5),
                                   n1234=ecount(D1_D2_D3_D4), n1235=ecount(D1_D2_D3_D5),
                                   n1245=ecount(D1_D2_D4_D5),  n1345=ecount(D1_D3_D4_D5),
                                   n2345=ecount(D2_D3_D4_D5),
                                   n12345=ecount(D1_D2_D3_D4_D5), category=c("D1", "D2", "D3", "D4", "D5"))

pdf(file="Network_overlap_for_D_replicates_150901_rahuls_network.pdf", width=6, height=6)
par(mfrow=c(1,2))
grid.draw(D_venn_edges)
grid.text(label="Zero-order Network Edges", x=0.8, y=0.9)
dev.off()

C1_C2 <- make_intersection(C1_zero, C2_zero)
C1_C3 <- make_intersection(C1_zero, C3_zero)
C1_C4 <- make_intersection(C1_zero, C4_zero)
C1_C5 <- make_intersection(C1_zero, C5_zero)
C2_C3 <- make_intersection(C2_zero, C3_zero)
C2_C4 <- make_intersection(C2_zero, C4_zero)
C2_C5 <- make_intersection(C2_zero, C5_zero)
C3_C4 <- make_intersection(C3_zero, C4_zero)
C3_C5 <- make_intersection(C3_zero, C5_zero)
C4_C5 <- make_intersection(C4_zero, C5_zero)
C1_C2_C3 <- make_intersection(C1_C2, C3_zero)
C1_C2_C4 <- make_intersection(C1_C2, C4_zero)
C1_C2_C5 <- make_intersection(C1_C2, C5_zero)
C1_C3_C4 <- make_intersection(C1_C3, C4_zero)
C1_C3_C5 <- make_intersection(C1_C3, C5_zero)
C1_C4_C5 <- make_intersection(C1_C4, C5_zero)
C2_C3_C4 <- make_intersection(C2_C3, C4_zero)
C2_C3_C5 <- make_intersection(C2_C3, C5_zero)
C2_C4_C5 <- make_intersection(C2_C4, C5_zero)
C3_C4_C5 <- make_intersection(C3_C4, C5_zero)
C1_C2_C3_C4 <- make_intersection(C1_C2, C3_C4)
C1_C2_C3_C5 <- make_intersection(C1_C2, C3_C5)
C1_C2_C4_C5 <- make_intersection(C1_C2, C4_C5)
C1_C3_C4_C5 <- make_intersection(C1_C3, C4_C5)
C2_C3_C4_C5 <- make_intersection(C2_C3, C4_C5)
C1_C2_C3_C4_C5 <- make_intersection(C1_C2_C3, C4_C5)


C_venn_edges <- draw.quintuple.venn(area1=ecount(C1_zero), area2=ecount(C2_zero), area3=ecount(C3_zero),
                                    area4=ecount(C4_zero), area5=ecount(C5_zero),
                                    n12=ecount(C1_C2), n13=ecount(C1_C3),
                                    n14=ecount(C1_C4), n15=ecount(C1_C5),
                                    n23=ecount(C2_C3), n24=ecount(C2_C4),
                                    n25=ecount(C2_C5), n34=ecount(C3_C4),
                                    n35=ecount(C3_C5), n45=ecount(C4_C5),
                                    n123=ecount(C1_C2_C3), n124=ecount(C1_C2_C4),
                                    n125=ecount(C1_C2_C5), n134=ecount(C1_C3_C4),
                                    n135=ecount(C1_C3_C5), n145=ecount(C1_C4_C5),
                                    n234=ecount(C2_C3_C4), n235=ecount(C2_C3_C5),
                                    n245=ecount(C2_C4_C5), n345=ecount(C3_C4_C5),
                                    n1234=ecount(C1_C2_C3_C4), n1235=ecount(C1_C2_C3_C5),
                                    n1245=ecount(C1_C2_C4_C5),  n1345=ecount(C1_C3_C4_C5),
                                    n2345=ecount(C2_C3_C4_C5),
                                    n12345=ecount(C1_C2_C3_C4_C5), category=c("C1", "C2", "C3", "C4", "C5"))

pdf(file="Network_overlap_for_C_replicates_150901_rahuls_network.pdf", width=6, height=6)
par(mfrow=c(1,2))
grid.draw(C_venn_edges)
grid.text(label="Zero-order Network Edges", x=0.8, y=0.9)
dev.off()



plot(D1_D2, vertex.label="", vertex.size=3)

D1_D2_vertex_vector <- which(V(D1_zero)$name%in%V(D1_D2)$name)
D1_D3_vertex_vector <- which(V(D1_zero)$name%in%V(D1_D3)$name)
D1_D4_vertex_vector <- which(V(D1_zero)$name%in%V(D1_D4)$name)
D1_D5_vertex_vector <- which(V(D1_zero)$name%in%V(D1_D5)$name)
plot(D1_zero, edge.color="grey", vertex.size=1, vertex.color="black",
     vertex.label="", mark.border=c("red", "blue", "green", "gold"), 
     mark.col=rgb(0, 0, 0, 0, maxColorValue=255), 
     mark.groups=list(D1_D2_vertex_vector, D1_D3_vertex_vector,
                      D1_D4_vertex_vector, D1_D5_vertex_vector))

#################################################
# Save lists of genes in each gene-list network #
# and genes in each overlap network             #
#################################################

setwd('/Users/pgriffin/Documents/Drosophila Selection Experiment/network_building_and_overlap_testing')

for(i in c('D1_zero', 'D2_zero', 'D3_zero', 'D4_zero', 'D5_zero',
           'C1_zero', 'C2_zero', 'C3_zero', 'C4_zero', 'C5_zero')){
  name_table <- data.frame(V(get(i))$name)
  file_name <- paste(i, 'order_network_gene_list.txt', sep="_")
  write.table(name_table, file=file_name, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

x <- c("D1", "D2", "D3", "D4", "D5")
combn2 <- combn(x, 2)
combn2 <- apply(combn2, 2, paste, collapse="_")
combn3 <- combn(x, 3)
combn3 <- apply(combn3, 2, paste, collapse="_")
combn4 <- combn(x, 4)
combn4 <- apply(combn4, 2, paste, collapse="_")


for(i in c('D1_D2', 'D2_zero', 'D3_zero', 'D4_zero', 'D5_zero')){
  name_table <- data.frame(V(get(i))$name)
  file_name <- paste(i, 'order_network_gene_list.txt', sep="_")
  write.table(name_table, file=file_name, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

for(i in c(combn2, combn3, combn4)){
  name_table <- data.frame(V(get(i))$name)
  if(nrow(name_table)>0){
    file_name <- paste(i, 'order_network_gene_list.txt', sep="_")
    write.table(name_table, file=file_name, row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}

###########################################
# Want to investigate for each network    #
# gene list: which of the genes belong    #
# to each body part / dev stage category? #
###########################################

setwd("~/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment/Gowinda_files")
bodypart_75pc <- read.table('Pmax_bodypart_75pc_filter.txt', sep="\t", stringsAsFactors=FALSE)
devstage_75pc <- read.table('Pmax_devstage_75pc_filter.txt', sep="\t", stringsAsFactors=FALSE)

x <- c("D1", "D2", "D3", "D4", "D5")
for(j in x){
 temp_network <- get(paste(j, "zero", sep="_"))  
 #set the layout here to stop it randomising positions for each plot
 testlayout <- layout.fruchterman.reingold(temp_network)
 file_name <- paste(j, "zero", "network_with_body_parts_highlighted.pdf", sep="_")
 pdf(file=file_name, width=8, height=8)
 for(i in 1:nrow(bodypart_75pc)){
    temp_category_list <- as.vector(strsplit(bodypart_75pc[i,3], split=" ")[[1]])
    temp_title <- paste(j, "zero", bodypart_75pc[i,1], sep="_")
    ### to plot as 'hairball'
    node_colours <- rep("black", times=length(V(temp_network)$name))
    edge_colours <- rgb(100, 100, 100, 20, maxColorValue=255)
    to_highlight <- which(V(temp_network)$name%in%temp_category_list)
    node_colours[to_highlight] <- "red"
    plot(temp_network, layout=testlayout, main=temp_title, 
         vertex.label="", vertex.color=node_colours, 
         vertex.size=2, add=FALSE)
    ### to plot as bipartite graph with highlighted genes separated
    #highlighted <- V(temp_network)$name%in%temp_category_list
    #plot(temp_network, layout=layout.bipartite(temp_network, types=highlighted), 
    #     main=temp_title, edge.color=edge_colours,
    #     vertex.label="", vertex.color=node_colours, 
    #     vertex.size=2, add=FALSE)
   }
 dev.off()
}

for(j in c(combn2, combn3, combn4)){
  temp_network <- get(j)
  if(length(V(temp_network)$name)>0){
    testlayout <- layout.fruchterman.reingold(temp_network)
    file_name <- paste(j, "network_with_body_parts_highlighted.pdf", sep="_")
    pdf(file=file_name, width=24, height=16)
    par(mfrow=c(4,6))
    for(i in 1:nrow(bodypart_75pc)){
      temp_category_list <- as.vector(strsplit(bodypart_75pc[i,3], split=" ")[[1]])
      temp_title <- paste(j, bodypart_75pc[i,1], sep="_")
      ### to plot as 'hairball'
      node_colours <- rep("black", times=length(V(temp_network)$name))
      edge_colours <- rgb(100, 100, 100, 20, maxColorValue=255)
      to_highlight <- which(V(temp_network)$name%in%temp_category_list)
      node_colours[to_highlight] <- "red"
      plot(temp_network, layout=testlayout, main=temp_title, 
           vertex.label="", vertex.color=node_colours, vertex.frame.color=node_colours,
           vertex.size=2, add=FALSE)
    }
    dev.off()
  }
}


for(j in c(combn2, combn3, combn4)){
  temp_network <- get(j)
  if(length(V(temp_network)$name)>0){
    testlayout <- layout.fruchterman.reingold(temp_network)
    file_name <- paste(j, "network_with_dev_stage_highlighted.pdf", sep="_")
    pdf(file=file_name, width=24, height=16)
    par(mfrow=c(5,6))
    for(i in 1:nrow(devstage_75pc)){
      temp_category_list <- as.vector(strsplit(devstage_75pc[i,3], split=" ")[[1]])
      temp_title <- paste(j, devstage_75pc[i,1], sep="_")
      ### to plot as 'hairball'
      node_colours <- rep("black", times=length(V(temp_network)$name))
      edge_colours <- rgb(100, 100, 100, 20, maxColorValue=255)
      to_highlight <- which(V(temp_network)$name%in%temp_category_list)
      node_colours[to_highlight] <- "red"
      plot(temp_network, layout=testlayout, main=temp_title, 
           vertex.label="", vertex.color=node_colours, vertex.frame.color=node_colours,
           vertex.size=4, add=FALSE)
    }
    dev.off()
  }
}

#########################################
# Resampling from full D. mel gene list #
#########################################

setwd("~/Documents/Drosophila Selection Experiment/genome_6.01")
gtf <- read.table('dmel-all-r6.01.gtf', sep="\t", stringsAsFactors=FALSE)
pre_gene_id <- sapply(strsplit(gtf[,9], split=";"), "[[", 1)
dup_gene_id <- sapply(strsplit(pre_gene_id, split=" "), "[[", 2)
gene_ids <- unique(dup_gene_id)

#make an array of 5 x 1000
# each element contains a resampled set of genes
# to appropriate size

to_resample <- list(D1_list, D2_list, D3_list, D4_list, D5_list)
names(to_resample) <- c("D1", "D2", "D3", "D4", "D5")

# resample_table=list()
# nIter=1000
# 
# for (i in 1:5) {
#   temp_real <- to_resample[[i]]
#   resample_table[[paste(names(to_resample[i]))]] = matrix(NA,nrow=length(temp_real))[,-1]
#   resample_table[[paste(names(to_resample[i]), 'networks', sep='_')]] = as.list(rep(NA, length=nIter))
#   for (j in 1:nIter) {
#     if(j%%100==0){
#       print(paste(names(to_resample[i]), " Resampling iteration #", j, sep=""))
#     }
#     temp_resample <- sample(gene_ids, size=length(temp_real), replace = FALSE)
#     resample_table[[paste(names(to_resample[i]))]]<-cbind(resample_table[[paste(names(to_resample[i]))]],temp_resample)
#     temp_graph <- make_zero_degree_subgraph(test_network, temp_resample)
#     resample_table[[paste(names(to_resample[i]), 'networks', sep='_')]][[j]] <- list(vcount=vcount(temp_graph), ecount=ecount(temp_graph), edgelist=get.edgelist(temp_graph),
#                                                                                      number_genes_in_background=length(which(temp_resample%in%V(test_network)$name)),
#                                                                                      density=graph.density(temp_graph, loops=FALSE))
#   }
# }
# 
# pdf("Zero-order network observed size, gene list length and density vs simulated 150707.pdf", width=16, height=16)
# par(mfcol=c(5,4))
# for(i in 1:5){
#   temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
#   temp_networks <- paste(temp_label, "_networks", sep="")
#   temp_observed_graph <- paste(temp_label, "_zero", sep="")
#   temp_vcount <- vcount(get(temp_observed_graph))
#   temp_dist <- sapply(get(temp_networks, resample_table), "[[", 1)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1] 
#   offset <- sd(temp_dist)*2.5
#   hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4), 
#        main=temp_label, xlab="Network size (no. nodes)",
#        ylim=c(0, 350))
#   lines(mydensity)
#   lines(x=c(temp_vcount, temp_vcount), y=c(0, 250), col="red")
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_vcount-50, y=100, labels=signif(normless, 4))
#   text(x=mean(temp_dist)-offset, y=325, labels=LETTERS[i], cex=2)
# }
# for(i in 1:5){
#   temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
#   temp_networks <- paste(temp_label, "_networks", sep="")
#   temp_observed_graph <- paste(temp_label, "_zero", sep="")
#   temp_ecount <- ecount(get(temp_observed_graph))
#   temp_dist <- sapply(get(temp_networks, resample_table), "[[", 2)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2.5
#   hist(temp_dist, main=temp_label, xlab="Network size (no. edges)",
#        ylim=c(0, 350), xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4))
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_ecount-50, y=100, labels=signif(normless,4))
#   lines(x=c(temp_ecount, temp_ecount), y=c(0, 200), col="red")
#   text(x=mean(temp_dist)-offset, y=325, labels=LETTERS[i+5], cex=2)
# }
# for(i in 1:5){
#   temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
#   temp_networks <- paste(temp_label, "_networks", sep="")
#   temp_observed_graph <- paste(temp_label, "_zero", sep="")
#   temp_number <- length(which(get(paste(temp_label, "_list", sep=""))%in%V(test_network)$name))
#   temp_dist <- sapply(get(temp_networks, resample_table), "[[", 4)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2.5
#   hist(temp_dist, main=temp_label, xlab="Number of genes with known interactions",
#        ylim=c(0, 400), xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*8))
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_number, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_number-50, y=100, labels=signif(normless,4))
#   lines(x=c(temp_number, temp_number), y=c(0, 200), col="red")
#   text(x=mean(temp_dist)-offset, y=375, labels=LETTERS[i+10], cex=2)
# }
# for(i in 1:5){
# ## This plots the graph density (average connectivity over all nodes)
# ## but the resampled gene lists aren't really the appropriate null
# ## comparison for this test. Instead, comparing to a network
# ## with the same number of nodes (see below...)
# ## - simulations found in resample_table_2
#   temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
#   temp_networks <- paste(temp_label, "_length_networks", sep="")
#   temp_observed_graph <- paste(temp_label, "_zero", sep="")
#   temp_density <- graph.density(get(temp_observed_graph), loops=FALSE)
#   temp_dist <- sapply(get(temp_networks, resample_table_2), "[[", 4)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2.5
#   hist(temp_dist, main=temp_label, xlab="Graph density (average connectivity over all nodes)",
#        ylim=c(0, 350), xlim=c(mean(temp_dist)-offset*2, mean(temp_dist)+offset*2))
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_density, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_density, y=250, labels=signif(normless,4))
#   lines(x=c(temp_density, temp_density), y=c(0, 200), col="red")
#   text(x=mean(temp_dist)-offset*2, y=325, labels=LETTERS[i+15], cex=2)
# }
# dev.off()

######################################
# Using resampled gene list networks #
# to calculate overlap               #
######################################

# resample_overlap_calcs=list()
# x <- c("D1", "D2", "D3", "D4", "D5")
# combn2 <- combn(x, 2)
# combn3 <- combn(x, 3)
# combn4 <- combn(x, 4)
# 
# for (i in 1:ncol(combn2)){
#   name1 <- combn2[1,i]
#   name2 <- combn2[2,i]
#   resample_overlap_calcs[[paste(name1, name2, sep="_")]]=list()
#   for (j in 1:nIter){
#     if(j%%100==0){
#       print(paste(name1, name2, "Overlap testing iteration #", j, sep=" "))
#     }
#     temp_graph1 <- graph.data.frame(resample_table[[paste(name1, 'networks', sep='_')]][[j]][[3]])
#     temp_graph2 <- graph.data.frame(resample_table[[paste(name2, 'networks', sep='_')]][[j]][[3]])
#     temp_overlap <- make_intersection(temp_graph1, temp_graph2)
#     resample_overlap_calcs[[paste(name1, name2, sep="_")]][[j]] <- list(vcount=vcount(temp_overlap), ecount=ecount(temp_overlap), edgelist=get.edgelist(temp_overlap))
#   }
# }
# for (i in 1:ncol(combn3)){
#   name1 <- combn3[1,i]
#   name2 <- combn3[2,i]
#   name3 <- combn3[3,i]
#   resample_overlap_calcs[[paste(name1, name2, name3, sep="_")]]=list()
#   for (j in 1:nIter){
#     if(j%%100==0){
#       print(paste(name1, name2, name3, "Overlap testing iteration #", j, sep=" "))
#     }
#     temp_graph1 <- graph.data.frame(resample_table[[paste(name1, 'networks', sep='_')]][[j]][[3]])
#     temp_graph2 <- graph.data.frame(resample_table[[paste(name2, 'networks', sep='_')]][[j]][[3]])
#     temp_graph3 <- graph.data.frame(resample_table[[paste(name3, 'networks', sep="_")]][[j]][[3]])
#     pre_temp_overlap <- make_intersection(temp_graph1, temp_graph2)
#     temp_overlap <- make_intersection(pre_temp_overlap, temp_graph3)
#     resample_overlap_calcs[[paste(name1, name2, name3, sep="_")]][[j]] <- list(vcount=vcount(temp_overlap), ecount=ecount(temp_overlap), edgelist=get.edgelist(temp_overlap))
#   }
# }
# 
# pdf("Zero-order network observed overlap vs simulated overlap 150703.pdf", width=20, height=32)
# par(mfcol=c(10,4))
# for(i in 1:ncol(combn2)){
#   name1 <- combn2[1,i]
#   name2 <- combn2[2,i]
#   temp_observed_overlap <- paste(name1, name2, sep="_")
#   temp_sim_overlap <- paste(name1, name2, sep="_")
#   temp_vcount <- vcount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs), "[[", 1)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2
#   hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*10), 
#        main=temp_observed_overlap, xlab="Overlap size (no. nodes)",
#        ylim=c(0, 900))
#   lines(mydensity)
#   lines(x=c(temp_vcount, temp_vcount), y=c(0, 300), col="red")
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_vcount-offset/2, y=100, labels=signif(normless, 4))
#   text(x=mean(temp_dist)-offset, y=850, labels=LETTERS[i], cex=2)
# }
# for(i in 1:ncol(combn2)){
#   name1 <- combn2[1,i]
#   name2 <- combn2[2,i]
#   temp_observed_overlap <- paste(name1, name2, sep="_")
#   temp_sim_overlap <- paste(name1, name2, sep="_")
#   temp_ecount <- ecount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs), "[[", 2)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2
#   hist(temp_dist, main=temp_observed_overlap, 
#        xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*10),
#        xlab="Overlap size (no. edges)",
#        ylim=c(0, 900))
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_ecount-offset/2, y=100, labels=signif(normless,4))
#   lines(x=c(temp_ecount, temp_ecount), y=c(0, 300), col="red")
#   text(x=mean(temp_dist)-offset, y=850, labels=LETTERS[i+10], cex=2)
# }
# extraletters <- c(LETTERS[21:26], paste("A", LETTERS, sep=""))
# for(i in 1:ncol(combn3)){
#   name1 <- combn3[1,i]
#   name2 <- combn3[2,i]
#   name3 <- combn3[3,i]
#   temp_observed_overlap <- paste(name1, name2, name3, sep="_")
#   temp_sim_overlap <- paste(name1, name2, name3, sep="_")
#   temp_vcount <- vcount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs), "[[", 1)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2
#   hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*15),
#        main=temp_observed_overlap, xlab="Overlap size (no. nodes)",
#        ylim=c(0, 1100))
#   lines(mydensity)
#   lines(x=c(temp_vcount, temp_vcount), y=c(0, 500), col="red")
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=mean(temp_dist)+offset*4, y=100, labels=signif(normless, 4))
#   text(x=mean(temp_dist)-offset, y=1000, labels=extraletters[i], cex=2)
# }
# 
# for(i in 1:ncol(combn3)){
#   name1 <- combn3[1,i]
#   name2 <- combn3[2,i]
#   name3 <- combn3[3, i]
#   temp_observed_overlap <- paste(name1, name2, name3, sep="_")
#   temp_sim_overlap <- paste(name1, name2, name3, sep="_")
#   temp_ecount <- ecount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs), "[[", 2)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   offset <- sd(temp_dist)*2
#   hist(temp_dist, main=temp_observed_overlap, 
#        xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*15),
#        xlab="Overlap size (no. edges)", ylim=c(0, 1100))
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=mean(temp_dist)+offset*4, y=100, labels=signif(normless,4))
#   lines(x=c(temp_ecount, temp_ecount), y=c(0, 500), col="red")
#   text(x=mean(temp_dist)-offset, y=1000, labels=extraletters[i+10], cex=2)
# }
# dev.off()


#########################################
# Repeat this with the gene lists that  #
# came from resampling at the SNP level #
# while retaining LD                    #
#########################################

ld_gene_lists <- readRDS(file="/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/gene_list_overlap_testing/Resampled_gene_lists_from_SNPs_with_LD_FBgn_format.rds")

#to_resample <- list(D1_list, D2_list, D3_list, D4_list, D5_list)
name_vector <- c("C1", "C2", "C3", "C4", "C5", 
                 "D1", "D2", "D3", "D4", "D5")
column_vector <- c(1:5, 12:16)

nIter=1000

ld_resample_table <- list()
for (i in 1:10) {
  ld_resample_table[[i]] <- lapply(ld_gene_lists, "[[", column_vector[i])
}

for(i in 1:10){
  ld_resample_table[[paste(name_vector[i], 'networks', sep='_')]] <- list()
  for (j in 1:nIter) {
    if(j%%100==0){
      print(paste(name_vector[i], " Resampling iteration #", j, sep=""))
    }
    temp_resample <- ld_resample_table[[i]][[j]]
    temp_graph <- make_zero_degree_subgraph(test_network, temp_resample)
    ld_resample_table[[paste(name_vector[i], 'networks', sep='_')]][[j]] <- list(vcount=vcount(temp_graph), ecount=ecount(temp_graph), edgelist=get.edgelist(temp_graph),
                                                                                     number_genes_in_background=length(which(temp_resample%in%V(test_network)$name)),
                                                                                     density=graph.density(temp_graph, loops=FALSE))
  }
}

saveRDS(ld_resample_table, file="/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/network_building_and_overlap_testing/Networks_from_resampled_SNPs_rahuls_network.rds")

#making overlap calcs
ld_resample_overlap_calcs=list()
x <- c("D1", "D2", "D3", "D4", "D5")
combn2 <- combn(x, 2)
combn3 <- combn(x, 3)
combn4 <- combn(x, 4)

for (i in 1:ncol(combn2)){
  name1 <- combn2[1,i]
  name2 <- combn2[2,i]
  ld_resample_overlap_calcs[[paste(name1, name2, sep="_")]]=list()
  for (j in 1:nIter){
    if(j%%100==0){
      print(paste(name1, name2, "Overlap testing iteration #", j, sep=" "))
    }
    temp_graph1 <- graph.data.frame(ld_resample_table[[paste(name1, 'networks', sep='_')]][[j]][[3]])
    temp_graph2 <- graph.data.frame(ld_resample_table[[paste(name2, 'networks', sep='_')]][[j]][[3]])
    temp_overlap <- make_intersection(temp_graph1, temp_graph2)
    ld_resample_overlap_calcs[[paste(name1, name2, sep="_")]][[j]] <- list(vcount=vcount(temp_overlap), ecount=ecount(temp_overlap), edgelist=get.edgelist(temp_overlap))
  }
}
for (i in 1:ncol(combn3)){
  name1 <- combn3[1,i]
  name2 <- combn3[2,i]
  name3 <- combn3[3,i]
  ld_resample_overlap_calcs[[paste(name1, name2, name3, sep="_")]]=list()
  for (j in 1:nIter){
    if(j%%100==0){
      print(paste(name1, name2, name3, "Overlap testing iteration #", j, sep=" "))
    }
    temp_graph1 <- graph.data.frame(ld_resample_table[[paste(name1, 'networks', sep='_')]][[j]][[3]])
    temp_graph2 <- graph.data.frame(ld_resample_table[[paste(name2, 'networks', sep='_')]][[j]][[3]])
    temp_graph3 <- graph.data.frame(ld_resample_table[[paste(name3, 'networks', sep="_")]][[j]][[3]])
    pre_temp_overlap <- make_intersection(temp_graph1, temp_graph2)
    temp_overlap <- make_intersection(pre_temp_overlap, temp_graph3)
    ld_resample_overlap_calcs[[paste(name1, name2, name3, sep="_")]][[j]] <- list(vcount=vcount(temp_overlap), ecount=ecount(temp_overlap), edgelist=get.edgelist(temp_overlap))
  }
}

saveRDS(ld_resample_overlap_calcs, file="Network_overlap_from_resampled_SNPs_rahul_network.rds")

#########################################
# using randomly chosen                 #
# networks of the same size as my obs   #
# ones instead, to investigate overlap #
#########################################

#Only used for plotting connectivity

sizes_to_resample <- list(D1_zero, D2_zero, D3_zero, D4_zero, D5_zero)
names(sizes_to_resample) <- c("D1", "D2", "D3", "D4", "D5")

resample_table_2=list()
nIter=1000

for (i in 1:5) {
  temp_real <- sizes_to_resample[[i]]
  resample_table_2[[paste(names(sizes_to_resample[i]))]] = matrix(NA,nrow=vcount(temp_real))[,-1]
  resample_table_2[[paste(names(sizes_to_resample[i]), 'length_networks', sep='_')]] = as.list(rep(NA, length=nIter))
  background_vertices <- get.data.frame(test_network, what='vertices')[,1]
  for (j in 1:nIter) {
    if(j%%100==0){
      print(paste(names(sizes_to_resample[i]), " Resampling iteration #", j, sep=""))
    }
    temp_resample <- sample(background_vertices, size=vcount(temp_real), replace = FALSE)
    resample_table_2[[paste(names(sizes_to_resample[i]))]]<-cbind(resample_table_2[[paste(names(sizes_to_resample[i]))]],temp_resample)
    temp_graph <- make_zero_degree_subgraph(test_network, temp_resample)
    resample_table_2[[paste(names(sizes_to_resample[i]), 'length_networks', sep='_')]][[j]] <- list(vcount=vcount(temp_graph), 
                                                                                                    ecount=ecount(temp_graph), 
                                                                                                    edgelist=get.edgelist(temp_graph), 
                                                                                                    density=graph.density(temp_graph, loops=FALSE))
  }
}

resample_overlap_calcs_2=list()
x <- c("D1", "D2", "D3", "D4", "D5")
combn2 <- combn(x, 2)
combn3 <- combn(x, 3)
combn4 <- combn(x, 4)

for (i in 1:ncol(combn2)){
  name1 <- combn2[1,i]
  name2 <- combn2[2,i]
  resample_overlap_calcs_2[[paste(name1, name2, sep="_")]]=list()
  for (j in 1:nIter){
    if(j%%100==0){
      print(paste(name1, name2, "Overlap testing iteration #", j, sep=" "))
    }
    temp_graph1 <- graph.data.frame(resample_table_2[[paste(name1, 'length_networks', sep='_')]][[j]][[3]])
    temp_graph2 <- graph.data.frame(resample_table_2[[paste(name2, 'length_networks', sep='_')]][[j]][[3]])
    temp_overlap <- make_intersection(temp_graph1, temp_graph2)
    resample_overlap_calcs_2[[paste(name1, name2, sep="_")]][[j]] <- list(vcount=vcount(temp_overlap), ecount=ecount(temp_overlap), edgelist=get.edgelist(temp_overlap))
  }
}

saveRDS(resample_table_2, file="Networks_from_interacting_genes_only_rahul_network.rds")
saveRDS(resample_overlap_calcs_2, file="Network_overlap_from_interacting_genes_only_rahul_network.rds")


pdf("Zero-order network observed size, gene list length and density vs simulated 150901 rahul network.pdf", width=16, height=16)
par(mfcol=c(5,4))
for(i in 1:5){
  temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
  temp_networks <- paste(temp_label, "_networks", sep="")
  temp_observed_graph <- paste(temp_label, "_zero", sep="")
  temp_vcount <- vcount(get(temp_observed_graph))
  temp_list <- ld_resample_table[[temp_networks]]
  temp_dist <- sapply(temp_list, "[[", 1)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  offset <- sd(temp_dist)*2.5
  hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*2), 
       main=temp_label, xlab="Network size (no. nodes)",
       ylim=c(0, 350))
  lines(mydensity)
  lines(x=c(temp_vcount, temp_vcount), y=c(0, 250), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_vcount, y=275, labels=signif(normless, 4))
  text(x=mean(temp_dist)-offset, y=325, labels=LETTERS[i], cex=2)
}
for(i in 1:5){
  temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
  temp_networks <- paste(temp_label, "_networks", sep="")
  temp_observed_graph <- paste(temp_label, "_zero", sep="")
  temp_ecount <- ecount(get(temp_observed_graph))
  temp_list <- ld_resample_table[[temp_networks]]
  temp_dist <- sapply(temp_list, "[[", 2)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2.5
  hist(temp_dist, main=temp_label, xlab="Network size (no. edges)",
       ylim=c(0, 350), xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*2))
  lines(mydensity)
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_ecount, y=225, labels=signif(normless,4))
  lines(x=c(temp_ecount, temp_ecount), y=c(0, 200), col="red")
  text(x=mean(temp_dist)-offset, y=325, labels=LETTERS[i+5], cex=2)
}
for(i in 1:5){
  temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
  temp_networks <- paste(temp_label, "_networks", sep="")
  temp_observed_graph <- paste(temp_label, "_zero", sep="")
  temp_number <- length(which(get(paste(temp_label, "_list", sep=""))%in%V(test_network)$name))
  temp_list <- ld_resample_table[[temp_networks]]
  temp_dist <- sapply(temp_list, "[[", 4)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2.5
  hist(temp_dist, main=temp_label, xlab="Number of genes with known interactions",
       ylim=c(0, 350), xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*3))
  lines(mydensity)
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_number, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_number, y=225, labels=signif(normless,4))
  lines(x=c(temp_number, temp_number), y=c(0, 200), col="red")
  text(x=mean(temp_dist)-offset, y=375, labels=LETTERS[i+10], cex=2)
}
for(i in 1:5){
  ## This plots the graph density (average connectivity over all nodes)
  ## but the resampled gene lists aren't really the appropriate null
  ## comparison for this test. Instead, comparing to a network
  ## with the same number of nodes (see below...)
  ## - simulations found in resample_table_2
  temp_label <- c("D1", "D2", "D3", "D4", "D5")[i]
  temp_networks <- paste(temp_label, "_length_networks", sep="")
  temp_observed_graph <- paste(temp_label, "_zero", sep="")
  temp_density <- graph.density(get(temp_observed_graph), loops=FALSE)
  temp_dist <- sapply(get(temp_networks, resample_table_2), "[[", 4)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2.5
  hist(temp_dist, main=temp_label, xlab="Graph density (average connectivity over all nodes)",
       ylim=c(0, 350), xlim=c(mean(temp_dist)-offset*2, mean(temp_dist)+offset*2))
  lines(mydensity)
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_density, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_density, y=250, labels=signif(normless,4))
  lines(x=c(temp_density, temp_density), y=c(0, 200), col="red")
  text(x=mean(temp_dist)-offset*2, y=325, labels=LETTERS[i+15], cex=2)
}
dev.off()




pdf("Zero-order network observed overlap vs overlap simulated from LD SNPs 150901 rahul network.pdf", width=20, height=32)
par(mfcol=c(10,4))
for(i in 1:ncol(combn2)){
  name1 <- combn2[1,i]
  name2 <- combn2[2,i]
  temp_observed_overlap <- paste(name1, name2, sep="_")
  temp_sim_overlap <- paste(name1, name2, sep="_")
  temp_vcount <- vcount(get(temp_observed_overlap))
  temp_dist <- sapply(get(temp_sim_overlap, ld_resample_overlap_calcs), "[[", 1)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2
  hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4), 
       main=temp_observed_overlap, xlab="Overlap size (no. nodes)",
       ylim=c(0, 900))
  lines(mydensity)
  lines(x=c(temp_vcount, temp_vcount), y=c(0, 300), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_vcount, y=350, labels=signif(normless, 4))
  text(x=mean(temp_dist)-offset, y=850, labels=LETTERS[i], cex=2)
}
for(i in 1:ncol(combn2)){
  name1 <- combn2[1,i]
  name2 <- combn2[2,i]
  temp_observed_overlap <- paste(name1, name2, sep="_")
  temp_sim_overlap <- paste(name1, name2, sep="_")
  temp_ecount <- ecount(get(temp_observed_overlap))
  temp_dist <- sapply(get(temp_sim_overlap, ld_resample_overlap_calcs), "[[", 2)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2
  hist(temp_dist, main=temp_observed_overlap, 
       xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4),
       xlab="Overlap size (no. edges)",
       ylim=c(0, 900))
  lines(mydensity)
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_ecount, y=350, labels=signif(normless,4))
  lines(x=c(temp_ecount, temp_ecount), y=c(0, 300), col="red")
  text(x=mean(temp_dist)-offset, y=850, labels=LETTERS[i+10], cex=2)
}
extraletters <- c(LETTERS[21:26], paste("A", LETTERS, sep=""))
for(i in 1:ncol(combn3)){
  name1 <- combn3[1,i]
  name2 <- combn3[2,i]
  name3 <- combn3[3,i]
  temp_observed_overlap <- paste(name1, name2, name3, sep="_")
  temp_sim_overlap <- paste(name1, name2, name3, sep="_")
  temp_vcount <- vcount(get(temp_observed_overlap))
  temp_dist <- sapply(get(temp_sim_overlap, ld_resample_overlap_calcs), "[[", 1)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2
  hist(temp_dist, xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4),
       main=temp_observed_overlap, xlab="Overlap size (no. nodes)",
       ylim=c(0, 1100))
  lines(mydensity)
  lines(x=c(temp_vcount, temp_vcount), y=c(0, 500), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_vcount, y=550, labels=signif(normless, 4))
  text(x=mean(temp_dist)-offset, y=1000, labels=extraletters[i], cex=2)
}

for(i in 1:ncol(combn3)){
  name1 <- combn3[1,i]
  name2 <- combn3[2,i]
  name3 <- combn3[3, i]
  temp_observed_overlap <- paste(name1, name2, name3, sep="_")
  temp_sim_overlap <- paste(name1, name2, name3, sep="_")
  temp_ecount <- ecount(get(temp_observed_overlap))
  temp_dist <- sapply(get(temp_sim_overlap, ld_resample_overlap_calcs), "[[", 2)
  print(shapiro.test(temp_dist))
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1]  
  offset <- sd(temp_dist)*2
  hist(temp_dist, main=temp_observed_overlap, 
       xlim=c(mean(temp_dist)-offset, mean(temp_dist)+offset*4),
       xlab="Overlap size (no. edges)", ylim=c(0, 1100))
  lines(mydensity)
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
  normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=temp_ecount, y=600, labels=signif(normless,4))
  lines(x=c(temp_ecount, temp_ecount), y=c(0, 500), col="red")
  text(x=mean(temp_dist)-offset, y=1000, labels=extraletters[i+10], cex=2)
}
dev.off()




# pdf("Zero-order network observed overlap vs simulated same-node-number overlap 150703.pdf", width=20, height=32)
# par(mfcol=c(10,2))
# for(i in 1:ncol(combn2)){
#   name1 <- combn2[1,i]
#   name2 <- combn2[2,i]
#   temp_observed_overlap <- paste(name1, name2, sep="_")
#   temp_sim_overlap <- paste(name1, name2, sep="_")
#   temp_vcount <- vcount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs_2), "[[", 1)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   hist(temp_dist, xlim=c(min(temp_dist), temp_vcount+100), main=temp_observed_overlap, xlab="Overlap size (no. nodes)")
#   lines(mydensity)
#   lines(x=c(temp_vcount, temp_vcount), y=c(0, 200), col="red")
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_vcount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_vcount-50, y=100, labels=signif(normless, 4))
# }
# for(i in 1:ncol(combn2)){
#   name1 <- combn2[1,i]
#   name2 <- combn2[2,i]
#   temp_observed_overlap <- paste(name1, name2, sep="_")
#   temp_sim_overlap <- paste(name1, name2, sep="_")
#   temp_ecount <- ecount(get(temp_observed_overlap))
#   temp_dist <- sapply(get(temp_sim_overlap, resample_overlap_calcs_2), "[[", 2)
#   print(shapiro.test(temp_dist))
#   
#   multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
#   mydensity <- density(temp_dist)
#   mydensity$y <- mydensity$y * multiplier[1]  
#   hist(temp_dist, main=temp_observed_overlap, 
#        xlim=c(min(temp_dist), max(temp_ecount)+100), 
#        xlab="Overlap size (no. edges)")
#   lines(mydensity)
#   myx <- seq(min(temp_dist), max(temp_dist), length = 100)
#   normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
#   lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
#   normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
#   #pval <- (1-normless)*2
#   text(x=temp_ecount-50, y=100, labels=signif(normless,4))
#   lines(x=c(temp_ecount, temp_ecount), y=c(0, 200), col="red")
# }
# dev.off()
# 
# 
# 




myx <- seq(min(temp_dist), max(temp_dist), length = 100)
normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
lines(myx, normal * multiplier[1], col = "blue", lwd = 2)
normless <- pnorm(q = temp_ecount, mean=mean(temp_dist), sd=sd(temp_dist))
pval <- (1-normless)*2
text(x=temp_ecount-50, y=200, labels=round(pval, 2))

pdf("Zero-order network observed size vs simulated size 150728.pdf", width=10, height=8)
for(i in c("D1", "D2", "D3", "D4", "D5")){
  temp_networks <- paste(i, "_networks", sep="")
  temp_observed_graph <- paste(i, "_zero", sep="")
  temp_vcount <- vcount(get(temp_observed_graph))
  temp_dist <- sapply(get(temp_networks, ld_resample_table), "[[", 1)
  hist(temp_dist, xlim=c(min(temp_dist)-100, max(temp_dist+200)), main=i, xlab="Network size (no. nodes)")
  lines(x=c(temp_vcount, temp_vcount), y=c(0, 200), col="red")
}
dev.off()

# Save observed networks (for investigating KEGG pathways etc)

setwd('/Users/pgriffin/Documents/Drosophila Selection Experiment/network_building_and_overlap_testing/rahuls_network_results')
for(i in Sample_code){
  temp_network <- get(paste(i, "_zero", sep=""))
  temp_filename <- paste(i, "_network_expanded_background.graphml", sep="")
  write.graph(temp_network, file=temp_filename, format="graphml")
}



#############################
# Redo network calculations #
# using gene lists produced #
# by simulation at the SNP  #
# level, retaining LD       #
#############################

#Converting from gene name to FBgn format is a complicating factor
#Trying the biomart package

# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# library("biomaRt")
# browseVignettes("biomaRt")
# 
# ensembl = useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
# filters = listFilters(ensembl)
# attributes = listAttributes(ensembl)
# #attribute 41: flybase_gene_id
# #attribute 44: flybasename_gene
# testid <- "drongo"
# testbm <- getBM(attributes=c("flybase_gene_id", "flybasename_gene"), 
#                 filters="flybasename_gene", 
#                 values=testid,
#                 mart=ensembl)
# 
# rearrange_iterations <- readRDS(file="Resampled_gene_lists_from_SNPs_with_LD.rds")
# 
# testlist <- rearrange_iterations[[1]][[1]]
# testbm <- getBM(attributes=c("flybase_gene_id", "flybasename_gene"), 
#                 filters="flybasename_gene", 
#                 values=testlist,
#                 mart=ensembl)
# testlist[testlist%in%testbm[,2]==FALSE]

######################################
# Now want to figure out how to plot #
# overlap in a useful way            #
######################################

library(HiveR)

set_axis_to_degree <- function(hpd_object, degree_interval_vector){
  #Function to change the plotting axes to plot by degree
  axis_2_min <- degree_interval_vector[1]
  axis_2_max <- degree_interval_vector[2]
  
  for(each_row in 1:nrow(hpd_object[[1]])){
    temp_radius <- hpd_object[[1]][each_row, 4]
    if(temp_radius>axis_2_min&temp_radius<=axis_2_max){
      hpd_object[[1]][each_row, 3] <- 2L
    }
    else if(temp_radius>axis_2_max){
      hpd_object[[1]][each_row, 3] <- 3L
    }
  }
  #hpd_object[[1]][3] <- as.integer(hpd_object[[1]][3])
  return(hpd_object)
}

get_edge_names_from_igraph <- function(igraph_object){
  edge_df <- get.data.frame(igraph_object, what = "edges")
  edge_name <- paste(edge_df$from, edge_df$to, sep="_")
  return(edge_name)
}

add_edge_names_to_hpd <- function(hpd_object){
  edge_name_vector <- c()
  for(i in 1:nrow(hpd_object[[2]])){
    temp_edge_row <- hpd_object[[2]][i,]
    name1 <- hpd_object[[1]][temp_edge_row[,1],2]
    name2 <- hpd_object[[1]][temp_edge_row[,2],2]
    edge_name <- paste(name1, name2, sep="_")
    edge_name_vector <- c(edge_name_vector, edge_name)
  }
  hpd_object[[2]]$lab <- edge_name_vector
  return(hpd_object)
}

plot.igraph(D_union, vertex.label="", vertex.size=0.25, vertex.c)

D1_edgenames <- get_edge_names_from_igraph(D1_zero)
D2_edgenames <- get_edge_names_from_igraph(D2_zero)
D3_edgenames <- get_edge_names_from_igraph(D3_zero)
D4_edgenames <- get_edge_names_from_igraph(D4_zero)
D5_edgenames <- get_edge_names_from_igraph(D5_zero)
D_union_edgenames <- get_edge_names_from_igraph(D_union)

D_union_edge_pres_abs_2 <- data.frame(D_union_edgenames%in%D1_edgenames,
                                    D_union_edgenames%in%D2_edgenames,
                                    D_union_edgenames%in%D3_edgenames,
                                    D_union_edgenames%in%D4_edgenames,
                                    D_union_edgenames%in%D5_edgenames)
D_union_edge_pres_abs_count_2 <- rowSums(D_union_edge_pres_abs_2)
colour_for_hairball <- rep("lightgrey", times=length(D_union_edge_pres_abs_count_2))
colour_for_hairball[which(D_union_edge_pres_abs_count_2==5)] <- "red"
colour_for_hairball[which(D_union_edge_pres_abs_count_2==4)] <- "orange"
colour_for_hairball[which(D_union_edge_pres_abs_count_2==3)] <- "yellow"
colour_for_hairball[which(D_union_edge_pres_abs_count_2==2)] <- "green"
colour_for_hairball[which(D_union_edge_pres_abs_count_2==1)] <- "blue"

l <- layout.fruchterman.reingold(g,niter=500,area=vcount(g)^2.3,repulserad=vcount(g)^2.8)
plot(g,layout=l)

jpeg(file="Test graph plotting.jpeg", width=20, height=20, res=300, units="cm")
plot.igraph(D_union, vertex.label="", vertex.size=0.2, 
            edge.color=colour_for_hairball, edge.curved=TRUE, 
            layout=layout.fruchterman.reingold(D_union, niter=10000))
dev.off()

D1_cc <- test_cc <- transitivity(D1_zero, type="localundirected", isolates="zero")
D1_n <- neighborhood.size(D1_zero, order=1)-1
D1_nn <- neighborhood.size(D1_zero, order=2)-D1_n-1
D1_branching <- D1_nn/D1_n
D1_zero <- set.vertex.attribute(D1_zero, name="Clustering coefficient", value=D1_cc)
D1_zero <- set.vertex.attribute(D1_zero, name="Branching", value=D1_branching)
D1_edge__df <- get.data.frame(D1_zero, what="edges")
D1_edge__df[,3] <- 1
D1_hive1 <- edge2HPD(edge_df=D1_edge__df, type="2D")
D1_hive2 <- mineHPD(D1_hive1, option="rad <- tot.edge.count")
D1_hive3 <- set_axis_to_degree(D1_hive2, c(2,10))
D1_unsorted <- data.frame(get.vertex.attribute(D1_zero, name="name"),
                          get.vertex.attribute(D1_zero, name="Clustering coefficient"),
                          get.vertex.attribute(D1_zero, name="Branching"))
D1_unsorted[,1] <- ordered(D1_unsorted[,1], levels=D1_hive3[[1]]$lab)
D1_sorted <- D1_unsorted[order(D1_unsorted[,1]),]


# set radius to clustering coefficient

D1_hive3[[1]]$radius <- D1_sorted[,2]
D1_hive3$axis.cols=c("lightgrey", "lightgrey", "lightgrey")
D1_hive3[[1]]$color <- rgb(100, 100, 100, 255, maxColorValue=255)
D1_hive3[[2]]$color <- rgb(20, 20, 20, 20, maxColorValue=255)
D1_hive3[[1]]$size <- 0.2
D1_hive3 <- add_edge_names_to_hpd(D1_hive3)
D1_hive4 <- mineHPD(D1_hive3, option = "remove zero edge")

pdf("D1_zeroorder_network.pdf", width=10, height=10)
plotHive(D1_hive4, method = "abs", bkgnd = "white", axLab.pos = 0.2, ch=0.25, 
         axLabs=c("deg 1-2", "deg 3-9", "deg > 10"),
         axLab.gpar=gpar(col="black", fontsize=10, lwd=1))
dev.off()

plot3dHive(D1_hive4, method="abs")


D_union_cc <- transitivity(D_union, type="localundirected", isolates="zero")
D_union_n <- neighborhood.size(D_union, order=1)-1
D_union_nn <- neighborhood.size(D_union, order=2)-D_union_n-1
D_union_branching <- D_union_nn/D_union_n
D_union_degree <- degree(D_union)
D_union <- set.vertex.attribute(D_union, name="Clustering coefficient", value=D_union_cc)
D_union <- set.vertex.attribute(D_union, name="Branching", value=D_union_branching)
D_union <- set.vertex.attribute(D_union, name="Degree", value=D_union_degree)
D_union_edge__df <- get.data.frame(D_union, what="edges")
D_union_edge__df[,3] <- 1
D_union_hive1 <- edge2HPD(edge_df=D_union_edge__df, type="2D")
D_union_hive2 <- mineHPD(D_union_hive1, option="rad <- tot.edge.count")

####### CHOOSE AXIS VARIABLE #######

#set axis to degree (connectivity)
#D_union_hive3 <- set_axis_to_degree(D_union_hive2, c(8,25))

#set axis to number of gene lists node appears in
node_appearance_D1 <- D_union_hive3[[1]]$lab%in%D1_list
node_appearance_D2 <- D_union_hive3[[1]]$lab%in%D2_list
node_appearance_D3 <- D_union_hive3[[1]]$lab%in%D3_list
node_appearance_D4 <- D_union_hive3[[1]]$lab%in%D4_list
node_appearance_D5 <- D_union_hive3[[1]]$lab%in%D5_list
node_appearance <- data.frame(node_appearance_D1, node_appearance_D2,
                              node_appearance_D3, node_appearance_D4,
                              node_appearance_D5)
node_appearance_number <- rowSums(node_appearance)
D_union_hive3[[1]]$axis <- 1L
D_union_hive3[[1]]$axis[which(node_appearance_number%in%c(3,4,5))] <- 3L
D_union_hive3[[1]]$axis[which(node_appearance_number==2)] <- 2L

D_union_unsorted <- data.frame(get.vertex.attribute(D_union, name="name"),
                          get.vertex.attribute(D_union, name="Clustering coefficient"),
                          get.vertex.attribute(D_union, name="Branching"),
                          get.vertex.attribute(D_union, name="Degree"))
D_union_unsorted[,1] <- ordered(D_union_unsorted[,1], levels=D_union_hive3[[1]]$lab)
D_union_sorted <- D_union_unsorted[order(D_union_unsorted[,1]),]

####### CHOOSE RADIUS VARIABLE ########

# set radius to clustering coefficient
D_union_hive3[[1]]$radius <- D_union_sorted[,2]

# set radius to degree
D_union_hive3[[1]]$radius <- D_union_sorted[,4]

# set radius to branching
D_union_hive3[[1]]$radius <- D_union_sorted[,3]

D_union_hive3$axis.cols=c("lightgrey", "lightgrey", "lightgrey")
D_union_hive3[[1]]$color <- rgb(100, 100, 100, 255, maxColorValue=255)
D_union_hive3[[2]]$color <- rgb(20, 20, 20, 20, maxColorValue=255)
D_union_hive3[[1]]$size <- 0.2
D_union_hive3 <- add_edge_names_to_hpd(D_union_hive3)
D_union_hive4 <- mineHPD(D_union_hive3, option = "remove zero edge")

D1_edgenames <- get_edge_names_from_igraph(D1_zero)
D2_edgenames <- get_edge_names_from_igraph(D2_zero)
D3_edgenames <- get_edge_names_from_igraph(D3_zero)
D4_edgenames <- get_edge_names_from_igraph(D4_zero)
D5_edgenames <- get_edge_names_from_igraph(D5_zero)

D_union_edge_pres_abs <- data.frame(D_union_hive4[[2]]$lab%in%D1_edgenames,
                                    D_union_hive4[[2]]$lab%in%D2_edgenames,
                                    D_union_hive4[[2]]$lab%in%D3_edgenames,
                                    D_union_hive4[[2]]$lab%in%D4_edgenames,
                                    D_union_hive4[[2]]$lab%in%D5_edgenames)
D_union_edge_pres_abs_count <- rowSums(D_union_edge_pres_abs)

D_union_hive4[[2]]$color[which(D_union_edge_pres_abs_count==5)] <- rgb(255, 0, 0, 100, maxColorValue=255)
D_union_hive4[[2]]$color[which(D_union_edge_pres_abs_count==4)] <- rgb(255, 128, 0, 100, maxColorValue=255)
D_union_hive4[[2]]$color[which(D_union_edge_pres_abs_count==3)] <- rgb(255, 255, 0, 100, maxColorValue=255)
D_union_hive4[[2]]$color[which(D_union_edge_pres_abs_count==2)] <- rgb(102,178,255, 100, maxColorValue=255)
D_union_hive4[[2]]$weight[which(D_union_edge_pres_abs_count==5)] <- 2
D_union_hive4[[2]]$weight[which(D_union_edge_pres_abs_count==4)] <- 2

pdf("D_replicates_union_network_axis_as_replication_level_radius_as_branching.pdf", width=10, height=10)
plotHive(D_union_hive4, method = "abs", bkgnd = "white", axLab.pos = 10, ch=10, 
         axLabs=c("in 1", "in 2", "in 3-5"),
         axLab.gpar=gpar(col="black", fontsize=10, lwd=1))
dev.off()

#Now make plots that highlight each replicate in turn
#plotting radius as degree here


D_union_hive4[[2]]$color <- rgb(20, 20, 20, 20, maxColorValue=255)
D_union_hive4[[2]]$color[which(D_union_edge_pres_abs[,2])] <- rgb(0, 255, 0, 100, maxColorValue=255)

# finally reorder edges so that the highlighted ones get drawn on top
## NB this doesn't work that well, reordering doesn't happen as claimed
## Edges are actually drawn in the order 1-2, 2-3, 3-1, 1-3, 3-2, 2-1, 1-1, 2-2, 3-3
# (first and second axis respectively)
edges <- D_union_hive4$edges
edgesGreen12 <- subset(edges, color == rgb(0, 255, 0, 100, maxColorValue=255) & axis)
edgesGrey <- subset(edges, color == rgb(20, 20, 20, 20, maxColorValue=255))
edges <- rbind(edgesGrey, edgesGreen)
D_union_hive4$edges <- edges

pdf("D_replicates_union_network_axis_as_replication_level_radius_as_degree_D2_highlighted.pdf", width=10, height=10)
plotHive(D_union_hive4, method = "abs", bkgnd = "white", axLab.pos = 10, ch=10, 
         axLabs=c("in 1", "in 2", "in 3-5"),
         axLab.gpar=gpar(col="black", fontsize=10, lwd=1))
dev.off()
###############

v1 <- D_union_hive4[[1]]$lab[D_union_hive4[[1]]$axis==3L]
v2 <- D_union_hive4[[1]]$lab[D_union_hive4[[1]]$radius>0.9]
which(v1%in%v2)
v1[v1%in%v2]


#########################################
# Investigating first-order networks:   #
# want to see if they are tractable for #
# looking at pathway enrichment         #
#########################################

make_first_degree_subgraph <- function(background_network, gene_list){
  background_edge_list <- get.edgelist(background_network)
  #gene_present <- which(background_edge_list[,1]%in%gene_list|background_edge_list[,2]%in%gene_list)
  #pre_first_subgraph <- subgraph.edges(background_network, eids=gene_present)
  pre_first_subgraph <- induced.subgraph(graph=background_network, 
                                         vids=unlist(neighborhood(graph=background_network,order=1, 
                                                                  nodes=which(V(background_network)$name%in%gene_list))))
  nonzero_degree_vertices <- degree(pre_first_subgraph)>0
  numbering_2 <- 1:length(V(pre_first_subgraph)$name)
  numbers_to_keep_2 <- numbering_2[nonzero_degree_vertices]
  first_subgraph <- induced.subgraph(pre_first_subgraph, numbers_to_keep_2)
  first_subgraph <- simplify(first_subgraph, remove.multiple=TRUE, remove.loops=TRUE)
  edgelist <- get.edgelist(first_subgraph)
  edgelist2 <- data.frame(paste(edgelist[,1], edgelist[,2], sep="_"), paste(edgelist[,2], edgelist[,1], sep="_"))
  first_subgraph <- set.edge.attribute(first_subgraph, "name1", value=as.character(edgelist2[,1]))
  first_subgraph <- set.edge.attribute(first_subgraph, "name2", value=as.character(edgelist2[,2]))
  
  return(first_subgraph)
}

D1_first <- make_first_degree_subgraph(test_network, D1_list)

#Okay this makes huge networks, the one for D3 is almost the full background network.
#Ignore/abandon...




#Want to set node and edge colour based on whether the node is in our actual gene list
# or not

#seed_list <- as.character(read.table("D1_unreduced_network_seed_proteins.txt", header=TRUE)$Uniprot)

# colour_nodes_by_vector <- function(hpd_object, vertex_df, gene_list, colour_member, colour_nonmember){
#   colour_vector <- c()
#   for(i in 1:nrow(vertex_df)){
#     temp_id <- vertex_df[i,10]
#     if(temp_id%in%gene_list){
#       colour_vector <- c(colour_vector, colour_member)
#     }
#     else {colour_vector <- c(colour_vector, colour_nonmember)}
#   }
#   hpd_object[[1]]$color <- colour_vector
#   return(hpd_object)
# }
# 
# colour_edges_by_vector <- function(hpd_object, vertex_df,
#                                    gene_list, colour_both_nodes, 
#                                    colour_one_node, colour_neither_node){
#   colour_vector <- c()
#   number_and_id <- data.frame(as.character(row.names(vertex_df)), vertex_df$id)
#   colnames(number_and_id) <- c("number", "id")
#   member_numbers <- number_and_id[number_and_id$id%in%gene_list,]$number
#   edge_start <- as.character(hpd_object[[2]][1][,1])
#   edge_end <- as.character(hpd_object[[2]][2][,1])
#   first_is_member <- edge_start%in%member_numbers
#   second_is_member <- edge_end%in%member_numbers
#   for(i in 1:nrow(hpd_object[[2]])){
#     temp_row <- hpd_object[[2]][i,1:2]
#     if(first_is_member[i]==TRUE & second_is_member[i]==TRUE){
#       colour_vector <- c(colour_vector, colour_both_nodes)
#     }
#     else if((first_is_member[i]==TRUE & second_is_member[i]==FALSE) | first_is_member[i]==FALSE & second_is_member[i]==TRUE){
#       colour_vector <- c(colour_vector, colour_one_node)
#     }
#     else {colour_vector <- c(colour_vector, colour_neither_node)}
#   }
#   hpd_object[[2]]$color <- colour_vector
#   return(hpd_object)
# }

#hive5 <- colour_nodes_by_vector(hive4, vertex_df, seed_list, rgb(100, 0, 0, 100, maxColorValue=255),
#                                rgb(100, 100, 100, 100, maxColorValue=255))
#hive5[[1]]$size <- 0.5
#hive5[[1]]$color <- rgb(100, 100, 100, 100, maxColorValue=255)
#hive5[[2]]$weight <- 0.5
#hive5[[2]]$color <- rgb(100, 100, 100, 100, maxColorValue=255)

#hive6 <- colour_edges_by_vector(hive5, vertex_df, seed_list, 
#                                rgb(255, 0, 0, 100, maxColorValue=255), 
#                                rgb(0, 0, 255, 100, maxColorValue=255), 
#                                rgb(100, 100, 100, 100, maxColorValue=255))



# method="rank" might be useful, but pretty dense with such a big network
# method="invert", action=-1 displays axes in reverse




