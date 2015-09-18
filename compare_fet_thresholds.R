setwd("~/Documents/Drosophila Selection Experiment/raw_data_processing")
library(ggplot2)

thresholds <- read.table('tabulating_fet_thresholds.txt', sep="\t", header=TRUE, stringsAsFactors=FALSE)

plot1 <- ggplot(data=thresholds, aes(x=min_count, y=no_snps_called,
                                     colour=sample)) +
  geom_point() +
  geom_line()

plot2 <- ggplot(data=thresholds, aes(x=min_coverage, y=no_snps_called,
                                     colour=sample)) +
  geom_point() +
  geom_line() +
  facet_wrap(~min_count) +
  theme(legend.position = "none")

pdf("SNPs_called_versus_min_coverage.pdf", height=6, width=8)
plot2
dev.off()