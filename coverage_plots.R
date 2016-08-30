setwd("~/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists")

Sample_codes <- c("C1", "C2", "C3", "C4", "C5",
                  "D1", "D2", "D3", "D4", "D5")

C1_test <- read.table("C1_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
C2_test <- read.table("C2_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
C3_test <- read.table("C3_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
C4_test <- read.table("C4_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
C5_test <- read.table("C5_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")

D1_test <- read.table("D1_noC_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
D2_test <- read.table("D2_noC_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
D3_test <- read.table("D3_noC_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
D4_test <- read.table("D4_noC_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")
D5_test <- read.table("D5_noC_sig_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")

random <- read.table("random_SNPs_coverage_in_all_pools.txt", header=TRUE, sep="\t")

MB_random_coverage <- random[,"MB_counts"]

pdf(file="Coverage_in_candidate_SNPs.pdf", width=9, height=24)
par(mfcol=c(5,2))
for(i in Sample_codes){
  temp_sample <- i
  temp_number <- which(Sample_codes==i)
  temp_col <- ifelse(which(Sample_codes==i)<6, "blue", "red")
  temp_counts_colname <- paste(i, "_counts", sep="")
  temp_candidate_coverage_df <- get(paste(i, "_test", sep=""))
  temp_focal_candidate_coverage <- temp_candidate_coverage_df[,temp_counts_colname]
  temp_MB_candidate_coverage <- temp_candidate_coverage_df[,"MB_counts"]
  temp_random_coverage <- random[,temp_counts_colname]
  plot(density(temp_MB_candidate_coverage),
       xlim=c(0, 300), ylim=c(0, 0.12),
       main=i, xlab="Coverage depth")
  lines(density(MB_random_coverage), col="grey", lty=2)
  
  lines(density(temp_focal_candidate_coverage), col=temp_col)
  lines(density(temp_random_coverage, na.rm = TRUE), col=temp_col, lty=2)
  
  tt <- t.test(temp_focal_candidate_coverage, temp_random_coverage, alternative="less")
  ttp <- signif(tt$p.value, 2)
  #text(x=150, y=0.05, labels=paste("P =", ttp))
  text(x=5, y=0.114, labels=LETTERS[temp_number], cex=3)
}
dev.off()

library(tidyr)
library(dplyr)
library(lme4)

random_rs <- random %>% gather(Sample, Coverage, MB_counts:D5_counts)
D1_rs <- D1_test %>% gather(Sample, Coverage, MB_counts:D5_counts)
D1_rs$Locus_type <- "D1_candidate"

random_nD1 <- random[sample(1:nrow(random), nrow(D1_test)),]
random_nD1_rs <- random_nD1 %>% gather(Sample, Coverage, MB_counts:D5_counts)
random_nD1_rs$Locus_type <- "random"

random_and_D1 <- rbind(random_nD1_rs, D1_rs)
random_and_D1$Sample <- as.factor(random_and_D1$Sample)
random_and_D1$Locus_type <- as.factor(random_and_D1$Locus_type)

random_D1_lm <- glmer(Locus_type~Coverage+Sample+(1|chr.pos), data=random_and_D1,
                 subset = which(random_and_D1$chr.chr=="2L"), family="binomial")


plot(random$D1_counts~random$D2_counts, 
     col=rgb(255, 0, 0, 100, maxColorValue=255), 
     pch=19, cex=0.2, xlim=c(0, 200), ylim=c(0, 200))


# fit normal curve
x <- mtcars$mpg
h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2) 
