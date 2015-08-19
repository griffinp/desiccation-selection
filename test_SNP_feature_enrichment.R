library(ggplot2)
library(grid)

setwd("~/Documents/Drosophila Selection Experiment/snpeff_SNP_feature_enrichment")

dat <- read.table('SNP_feature_enrichment_table_for_R.txt', sep="\t", header=TRUE,
                  stringsAsFactors=FALSE)
for(i in 3:6){
  dat[,i] <- as.numeric(dat[,i])
}

dat$prop_all <- dat[,3]/dat[,4]
dat$prop_in <- dat[,5]/dat[,6]
dat$allout <- dat[,4]-dat[,3]
dat$sigout <- dat[,6]-dat[,5]

est <- c()
pval <- c()

for(i in 1:nrow(dat)){
  temp_row <- dat[i,]
  
  temp_mat <- matrix(c(temp_row$All.C1.in.category, temp_row$allout,
                       temp_row$Sig.C1.in.category, temp_row$sigout),
                     nrow=2, ncol=2)
  temp_test <- chisq.test(temp_mat)
  
  est <- c(est, temp_test$statistic)
  pval <- c(pval, temp_test$p.value)
}
results_table <- cbind(dat, est, pval)
write.table(results_table, file="SNP_feature_enrichment_chisq_results.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

#manually categorised categories into groups

chi <- read.table("SNP_feature_enrichment_chisq_results.txt", sep="\t",
                  stringsAsFactors=FALSE, header=TRUE)
chi <- subset(chi, Category%in%c('initiator_codon_variant ', 
                             'initiator_codon_variant+non_canonical_start_codon ',
                             'intragenic_variant ', 
                             'nc_transcript_variant ',
                             'start_lost ',
                             'stop_lost ',
                             'stop_retained ',
                             'NONE ',
                             'TRANSCRIPT ')==FALSE)


#for i in category groups
#for j in category

pdf(file="Plotting_SNP_categories.pdf", width=9, height=9)

for(i in 1:4){
  cat_group_level <- levels(as.factor(chi$Comparison_group))[i]
  temp_category_group <- chi[chi$Comparison_group==cat_group_level,]
  temp_category_group$pbin <- ifelse(temp_category_group$pval < 0.05, "black", "grey")
  all <- temp_category_group[,c(1,2,7,11,12)]
  sig <- temp_category_group[,c(1,2,8,11,12)]
  colnames(all)[3] <- "prop"
  colnames(sig)[3] <- "prop"
  
  temp <- rbind(all, sig)
  temp$from <- c(rep(c("all", "sig"), each=nrow(temp_category_group)))
  tempsig <- temp[temp$pval < 0.05,]
  
#   tempplot <- ggplot(data=temp, aes(x=Replicate, y=prop)) +
#     geom_point(data=tempsig, aes(x=Replicate, y=prop,
#                                  fill=Category, size=1),
#                #position=position_dodge(height=0, width=0.15), 
#                stat="identity") + 
#     geom_point(aes(shape=from, size=2, colour=Category), 
#                #position=position_dodge(height=0, width=0.15), 
#                stat="identity") +
#     #scale_colour_manual(values=c("red", "blue", "lightgreen")) + 
#     scale_shape(solid=FALSE, guide=FALSE) +
#     scale_size(guide=FALSE) +
#     #scale_shape(solid=TRUE) + 
#     theme_bw() +
#     theme(plot.background = element_blank(),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank()) +
#     theme(legend.title=element_blank(),
#           legend.position="none") +
#     facet_wrap(~Category)
  
  tempplot2 <- ggplot(data=temp_category_group) +
    geom_segment(aes(x=Replicate, y=prop_all, xend=Replicate, yend=prop_in,
                     colour=pbin), 
                 arrow = arrow(length=unit(0.1, "cm")),
                 stat="identity") + 
    scale_colour_manual(values=c("black", "grey")) + 
    #scale_shape(solid=FALSE, guide=FALSE) +
    #scale_size(guide=FALSE) +
    #scale_shape(solid=TRUE) + 
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    theme(legend.title=element_blank(),
          legend.position="none") +
    facet_wrap(~Category)
  
  
  plot(tempplot2)
}
dev.off()


i <- 2

cat_group_level <- levels(as.factor(chi$Comparison_group))[i]
temp_category_group <- chi[chi$Comparison_group==cat_group_level,]
temp_category_group$pbin <- ifelse(temp_category_group$pval < 0.05, "black", "grey")
temp_category_group$Category <- ordered(temp_category_group$Category,
                                        levels=c("LOW ", "MODERATE ", "HIGH ", "MODIFIER "))

tempplot2 <- ggplot(data=temp_category_group) +
  geom_segment(aes(x=Replicate, y=prop_all, xend=Replicate, yend=prop_in,
                   colour=pbin), 
               arrow = arrow(length=unit(0.1, "cm")),
               stat="identity") + 
  scale_colour_manual(values=c("black", "grey")) + 
  ylab("Change in proportion of features") +
  xlab("Replicate line") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank(),
        legend.position="none") +
  facet_wrap(~Category)

pdf(file="SNP_effect_size_enrichment.pdf", width=5, height=5)
plot(tempplot2)
dev.off()


i <- 3

cat_group_level <- levels(as.factor(chi$Comparison_group))[i]
temp_category_group <- chi[chi$Comparison_group==cat_group_level,]
temp_category_group$pbin <- ifelse(temp_category_group$pval < 0.05, "black", "grey")
#temp_category_group$Category <- ordered(temp_category_group$Category,
#                                        levels=c("LOW ", "MODERATE ", "HIGH ", "MODIFIER "))

tempplot4 <- ggplot(data=temp_category_group) +
  geom_segment(aes(x=Replicate, y=prop_all, xend=Replicate, yend=prop_in,
                   colour=pbin), 
               arrow = arrow(length=unit(0.1, "cm")),
               stat="identity") + 
  scale_colour_manual(values=c("black", "grey")) + 
  ylab("Change in proportion of features") +
  xlab("Replicate line") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank(),
        legend.position="none") +
  facet_wrap(~Category, ncol=2, nrow=10)

pdf(file="SNP_effect_location_enrichment.pdf", width=5, height=10)
plot(tempplot4)
dev.off()








